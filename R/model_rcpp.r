# SEIRS model with Rcpp implementation ----

## helper data and functions ----
### date data frame for assistance/reference
dates <- data.frame(date = seq(as.Date("23-03-2020", format = "%d-%m-%Y"), as.Date("02-03-2022", format = "%d-%m-%Y"), 1)) %>% 
  mutate(time = 0:(n()-1),
         fortnight = paste(isoyear(date), "/", sprintf("%02d", ceiling(isoweek(date)/2))),
         mmyyyy = format(date, "%m/%Y"),
         quarter = quarters(date)) %>% 
  mutate(fortnight_n = as.integer(factor(fortnight)))

### define function to choose which fortnightly contact matrix to use 
fortnight_lookup <- dates$fortnight_n

c_t <- function(t) {
  fortnight_matrix[[fortnight_lookup[findInterval(t, dates$time)]]]$matrix
}

# pre-flatten fortnight contact matrices into a C++-ready format.
# run after fortnight_matrix and fortnight_lookup created
#
# fortnight_matrix : list of length n_fn; each element has $matrix (9x9)
# fortnight_lookup : dates$fortnight_n - integer vector of length 710 mapping
#                   each day (1-indexed position = day 0-indexed + 1) to a
#                   fortnight index (1-based)
prepare_contacts_cpp <- function(fortnight_matrix, fortnight_lookup) {
  n_fn  <- length(fortnight_matrix)
  n_age <- 9L
  
  # Stack matrices row-major (transpose before as.vector because R is column-major)
  contacts_flat <- numeric(n_fn * n_age * n_age)
  for (i in seq_len(n_fn)) {
    mat <- fortnight_matrix[[i]]$matrix
    contacts_flat[((i - 1L) * 81L + 1L):(i * 81L)] <- as.vector(t(mat))
  }
  
  # fn_at_day[d] (0-indexed d = day) = fortnight index (1-based)
  # fortnight_lookup[findInterval(t, 0:709)] = fortnight_lookup[t+1] for t in 0:708
  # so fn_at_day is just fortnight_lookup as an integer vector
  fn_at_day <- as.integer(fortnight_lookup)
  
  list(contacts_flat = contacts_flat, fn_at_day = fn_at_day)
}

### create data frame defining model parameters
source(here("R", "create_combinations.r"))
combinations <- create_combinations()

### average population estimate between 2020-2022 for each age group
for (i in seq_len(nrow(scot_population))) {
  assign(paste0("tot", i), scot_population$average_n[i])
}

### define daily births and death rate
daily_births <- scot_births %>% select(births_daily) %>% pull() %>% mean() %>% floor() # load scot_births from explore_scotland.R
hazard_death <- daily_births/tot9 # for now keep death rate equal to birth rate, calculate per-capita death hazard

### define aging rates
a1 <- 1/(5*365) # 0-4
a2 <- 1/(6*365) # 5-10
a3 <- 1/(4*365) # 11-14
a4 <- 1/(10*365) # 15-24
a5 <- 1/(10*365) # 25-34
a6 <- 1/(10*365) # 35-44
a7 <- 1/(10*365) # 45-54
a8 <- 1/(10*365) # 55-64

# aging rate vector
a_vec <- c(a1, a2, a3, a4, a5, a6, a7, a8)

# C++ ODE implementation ----
sourceCpp(code = '
#include <Rcpp.h>
using namespace Rcpp;

// Compute SEIRS ODE right-hand side.
// State layout (interleaved): S1,E1,I1,R1, S2,E2,I2,R2, ..., S9,E9,I9,R9
// Contact matrix C is 9x9, row-major: C[i,j] = C_flat[i*9 + j]
static void compute_rhs(
    const double* y, double* dydt, double* inc_out,
    const double* C,
    double sigma, double gam, double omega,
    double p_inf, double bg, double mu_b, double mu_d,
    const double* a_vec, int n) {

  double S[9], E[9], I_[9], R_[9], N_[9];
  for (int i = 0; i < n; ++i) {
    S[i]  = y[4*i];
    E[i]  = y[4*i + 1];
    I_[i] = y[4*i + 2];
    R_[i] = y[4*i + 3];
    N_[i] = S[i] + E[i] + I_[i] + R_[i];
  }

  // Force of infection: lambda[i] = (C %*% (I/N))[i] * p_inf
  double lambda[9];
  for (int i = 0; i < n; ++i) {
    lambda[i] = 0.0;
    for (int j = 0; j < n; ++j) {
      double frac = (N_[j] > 0.0) ? I_[j] / N_[j] : 0.0;
      lambda[i] += C[i * n + j] * frac;
    }
    lambda[i] = lambda[i] * p_inf + bg;
  }

  // Derivatives
  for (int i = 0; i < n; ++i) {
    double a_out = (i < n - 1) ? a_vec[i] : mu_d;
    double S_in  = (i == 0) ? mu_b          : a_vec[i-1] * S[i-1];
    double E_in  = (i == 0) ? 0.0           : a_vec[i-1] * E[i-1];
    double I_in  = (i == 0) ? 0.0           : a_vec[i-1] * I_[i-1];
    double R_in  = (i == 0) ? 0.0           : a_vec[i-1] * R_[i-1];

    dydt[4*i]     = S_in - lambda[i]*S[i]  + omega*R_[i] - a_out*S[i];
    dydt[4*i + 1] = E_in + lambda[i]*S[i]  - sigma*E[i]  - a_out*E[i];
    dydt[4*i + 2] = I_in + sigma*E[i]      - gam*I_[i]   - a_out*I_[i];
    dydt[4*i + 3] = R_in + gam*I_[i]       - omega*R_[i] - a_out*R_[i];

    inc_out[i] = sigma * E[i];
  }
}

// RK4 ODE integration of the SEIRS model with time-varying contact matrices.
//
// Parameters:
//   y0            - initial state vector (length 36: S1,E1,I1,R1,...,S9,E9,I9,R9)
//   times         - integer day sequence to record output at (e.g. 0:709)
//   contacts_flat - all contact matrices stacked row-major (n_fn * 81)
//   fn_at_day     - fortnight index (1-based) for each day; fn_at_day[d] for day d (0-indexed)
//   sigma, gam    - incubation rate (1/inc_period), recovery rate (1/inf_period)
//   omega         - waning rate (1/imm_period)
//   p_inf         - probability of infection per contact
//   mu_b, mu_d    - daily birth flux and per-capita death hazard
//   a_vec         - aging rates for groups 1-8 (length 8)
//
// Returns a matrix with columns: time, S1,E1,I1,R1,...,S9,E9,I9,R9, inc1,...,inc9

// [[Rcpp::export]]
NumericMatrix seirs_ode_rcpp(
    NumericVector y0,
    IntegerVector times,
    NumericVector contacts_flat,
    IntegerVector fn_at_day,
    double sigma, double gam, double omega,
    double p_inf, double bg, double mu_b, double mu_d,
    NumericVector a_vec) {

  const int n_age    = 9;
  const int n_states = 36;
  const int n_cols   = 1 + n_states + n_age;  // time + states + incidences
  int n_times = times.size();
  int n_fn    = contacts_flat.size() / (n_age * n_age);
  int max_day = fn_at_day.size() - 1;

  NumericMatrix out(n_times, n_cols);

  // Working state vector
  std::vector<double> y(n_states);
  for (int i = 0; i < n_states; ++i) y[i] = y0[i];

  // RK4 temporary buffers
  std::vector<double> k1(n_states), k2(n_states), k3(n_states), k4(n_states);
  std::vector<double> ytmp(n_states), inc(n_age);

  // Retrieve pointer to contact matrix for integer day d
  auto get_C = [&](int d) -> const double* {
    if (d < 0)        d = 0;
    if (d > max_day)  d = max_day;
    int fn = fn_at_day[d] - 1;  // convert to 0-indexed
    if (fn < 0)       fn = 0;
    if (fn >= n_fn)   fn = n_fn - 1;
    return contacts_flat.begin() + fn * n_age * n_age;
  };

  for (int ti = 0; ti < n_times; ++ti) {
    int t = times[ti];

    // Record state and incidence at current time
    const double* C = get_C(t);
    compute_rhs(y.data(), k1.data(), inc.data(), C,
                sigma, gam, omega, p_inf, bg, mu_b, mu_d,
                a_vec.begin(), n_age);

    out(ti, 0) = t;
    for (int i = 0; i < n_states; ++i) out(ti, 1 + i)              = y[i];
    for (int i = 0; i < n_age;    ++i) out(ti, 1 + n_states + i)   = inc[i];

    // Integrate from t to times[ti+1] using daily RK4 steps
    if (ti < n_times - 1) {
      int n_steps = times[ti + 1] - t;

      for (int step = 0; step < n_steps; ++step) {
        int t_curr = t + step;
        int t_next = t_curr + 1;

        // k1 at t_curr
        C = get_C(t_curr);
        compute_rhs(y.data(), k1.data(), inc.data(), C,
                    sigma, gam, omega, p_inf, bg, mu_b, mu_d,
                    a_vec.begin(), n_age);

        // k2 at t_curr (same contact matrix; fortnights >> 1 day)
        for (int i = 0; i < n_states; ++i) ytmp[i] = y[i] + 0.5 * k1[i];
        compute_rhs(ytmp.data(), k2.data(), inc.data(), C,
                    sigma, gam, omega, p_inf, bg, mu_b, mu_d,
                    a_vec.begin(), n_age);

        // k3 at t_curr
        for (int i = 0; i < n_states; ++i) ytmp[i] = y[i] + 0.5 * k2[i];
        compute_rhs(ytmp.data(), k3.data(), inc.data(), C,
                    sigma, gam, omega, p_inf, bg, mu_b, mu_d,
                    a_vec.begin(), n_age);

        // k4 at t_next (may use a different fortnight contact matrix)
        const double* C_next = get_C(t_next);
        for (int i = 0; i < n_states; ++i) ytmp[i] = y[i] + k3[i];
        compute_rhs(ytmp.data(), k4.data(), inc.data(), C_next,
                    sigma, gam, omega, p_inf, bg, mu_b, mu_d,
                    a_vec.begin(), n_age);

        // Combine
        for (int i = 0; i < n_states; ++i) {
          y[i] += (1.0 / 6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
          if (y[i] < 0.0) y[i] = 0.0;  // numerical floor
        }
      }
    }
  }

  return out;
}
')

# run C++ SEIRS ODE solver.
# returns a matrix with the same column layout as deSolve::ode():
seirs_rcpp <- function(y0, times, parms, contacts_prepped, a_vec_in = a_vec) {
  out <- seirs_ode_rcpp(
    y0            = as.numeric(y0),
    times         = as.integer(times),
    contacts_flat = contacts_prepped$contacts_flat,
    fn_at_day     = contacts_prepped$fn_at_day,
    sigma         = parms$sigma,
    gam           = parms$gamma,
    omega         = parms$omega,
    p_inf         = parms$p_inf,
    bg           = parms$bg,
    mu_b          = parms$mu_b,
    mu_d          = parms$mu_d,
    a_vec         = a_vec_in
  )

  state_names <- paste0(rep(c("S", "E", "I", "R"), 9), rep(1:9, each = 4))
  colnames(out) <- c("time", state_names, paste0("inc", 1:9))

  out
}

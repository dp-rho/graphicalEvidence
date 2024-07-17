
/* Class for testing repeated looping code in isolation */

class FunctionTimer {
public:
  FunctionTimer(const std::string& timerName)
    : m_timerName(timerName)
  {}

  void TimerStart() {
    m_start = std::chrono::high_resolution_clock::now();
  }

  void TimerEnd() {
    auto end = std::chrono::high_resolution_clock::now();
    m_totalDuration += std::chrono::duration_cast<std::chrono::microseconds>(end - m_start).count();
  }

  void getTotalDuration(const int nruns) {
    long long divided_ms = m_totalDuration / nruns;
    std::cout << "Compiled execution time " << m_timerName << " (microseconds): " << divided_ms << std::endl;
  }

  void resetDuration() {
    m_totalDuration = 0;
  }

private:
  std::string m_timerName;
  std::chrono::time_point<std::chrono::high_resolution_clock> m_start;
  long long m_totalDuration;
};

class FunctionTimerArr {
public:
  FunctionTimerArr(const std::string& timerName, const unsigned int size)
    : m_timerName(timerName), m_size(size)
  {}

  void TimerStart() {
    m_start = std::chrono::high_resolution_clock::now();
  }

  void TimerEnd(const unsigned int index) {
    auto end = std::chrono::high_resolution_clock::now();
    m_totalDuration[index] += std::chrono::duration_cast<std::chrono::microseconds>(end - m_start).count();
  }

  void getTotalDuration(const int nruns) {
    for (unsigned int i = 0; i < m_size; i++) {
      long long divided_ms = m_totalDuration[i] / nruns;
      std::cout << "Compiled execution time of " << m_timerName << " " << i << " (microseconds): " << divided_ms << std::endl;
    }
  }

  void resetDuration() {
    for (unsigned int i = 0; i < m_size; i++) {
      m_totalDuration[i] = 0;
    }
  }

private:
  std::string m_timerName;
  const unsigned int m_size;
  std::chrono::time_point<std::chrono::high_resolution_clock> m_start;
  long long m_totalDuration[50];
};

/* Global timer variables */
extern FunctionTimer g_eq_9_timer;
extern FunctionTimer g_mcmc_hw_timer;
extern FunctionTimer g_eq_11_timer;
extern FunctionTimer g_mcmc_last_col_timer;
extern FunctionTimer g_inv_c_hw;
extern FunctionTimer g_mu_reduced1_hw;
extern FunctionTimer g_mu_reduced2_hw;
extern FunctionTimer g_mu_reduced3_hw;
extern FunctionTimer g_update_omega_hw1;
extern FunctionTimer g_update_omega_hw2;
extern FunctionTimer g_sample_omega_hw;
extern FunctionTimer g_inv_omega_11_hw;
extern FunctionTimer g_last_col_t1;
extern FunctionTimer g_last_col_t2;
extern FunctionTimer g_last_col_t3;
extern FunctionTimer g_last_col_t4;
extern FunctionTimer g_last_col_t5;
extern FunctionTimer g_last_col_t6;
extern FunctionTimer g_last_col_t7;
extern FunctionTimer g_last_col_t8;
extern FunctionTimer g_last_col_t9;

extern FunctionTimerArr g_last_col_solve;

extern double last_col_conds[50];
extern double last_col_min_conds[50];
extern double max_calc_time;
extern arma::mat inv_c_max_time;
extern arma::vec solve_for_max_time;

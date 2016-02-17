#ifndef _TMODPARAMETERS_HH_
#define _TMODPARAMETERS_HH_

extern "C" {
  void __modparameters_MOD_sethiggsmasswidth(double *mass, double *width);
  void __modparameters_MOD_setdecaymodes(int idfirst[2], int idsecond[2]);

  // NOTE: LOGICAL==INT, LOGICAL*1==BOOL(!???): http://www.yolinux.com/TUTORIALS/LinuxTutorialMixingFortranAndC.html
  void __modparameters_MOD_setspinzerovvcouplings        (double vvcoupl[39][2], int cqsq[3], double Lambda_qsq[4][3], int* usewwcoupl); // YES, THE LAST ARGUMENT IS AN INT!
  void __modparameters_MOD_setspinzerovvcouplings_nogamma(double vvcoupl[32][2], int cqsq[3], double Lambda_qsq[4][3], int* usewwcoupl); // YES, THE LAST ARGUMENT IS AN INT!
  void __modparameters_MOD_setdistinguishwwcouplingsflag(int* doallow); // YES, THE ARGUMENT IS AN INT!
  void __modparameters_MOD_setspinzeroggcouplings(double ggcoupl[3][2]);
  void __modparameters_MOD_setspinzeroqqcouplings(double qqcoupl[2][2]);

  void __modparameters_MOD_setspinonecouplings(double qqcoupl[2][2], double vvcoupl[2][2]);
  void __modparameters_MOD_setspintwocouplings(double acoupl[5][2], double bcoupl[10][2], double qlr[2][2]);
}

#endif


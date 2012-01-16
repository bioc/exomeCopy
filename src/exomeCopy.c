#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void negloglike( int *ptmax, int *pnstates, double *start_probs, double *a, double *emit_probs, double *alpha, double *alpha_new, double *nll )
{
  int tmax = ptmax[0];
  int nstates = pnstates[0];
  int i, j, t;
  double c = 0;
  double c_new, c_min, trans_probs_sum;
  
  R_CheckUserInterrupt();

  for (i=0; i<nstates; i++) {
    alpha[i] = start_probs[i] * emit_probs[i];
    c += alpha[i];
  }
  for (i=0; i<nstates; i++) {
    alpha[i] = alpha[i]/c;
  }
  c_min = c;
  c = log(c);
  for (t=1; t<tmax; t++) {
    R_CheckUserInterrupt();
    c_new = 0;
    for (i=0; i<nstates; i++) {
      trans_probs_sum = 0;
      for (j=0; j<nstates; j++) {
	trans_probs_sum += alpha[j]*a[j+i*nstates];
      }
      alpha_new[i] = emit_probs[i+t*nstates] * trans_probs_sum;
      c_new += alpha_new[i];
    }
    if (c_new > 0) {
      for (i=0; i<nstates; i++) {
	alpha[i] = alpha_new[i]/c_new;
      }
      c_min = fmin2(c_new,c_min);
    } else {
      /* assign the lowest seen likelihood */
      c_new = c_min;
    }
    c += log(c_new);
  }

  nll[0] = -1.0 * c;
}

void viterbi( int *ptmax, int *pnstates, double *start_probs, double *a, double *emit_probs, double *v, int *v_path, int *path, double *trans_probs, double *trans_probs_max, int *trans_probs_whichmax)
{
  int i, j, t;
  int tmax = ptmax[0];
  int nstates = pnstates[0];
  double vsum = 0.0;
  double running_max = 0.0;

  R_CheckUserInterrupt();

  for (i=0; i<nstates; i++) {
    v[i] = start_probs[i] * emit_probs[i];
  }
  for (t=1; t<tmax; t++) {
    R_CheckUserInterrupt();
    for (j=0; j<nstates; j++) {
      for (i=0; i<nstates; i++) {
	trans_probs[i+nstates*j] = a[i+nstates*j] * v[i+nstates*(t-1)];
	if (i==0) {
	  trans_probs_max[j] = trans_probs[i+nstates*j];
	  trans_probs_whichmax[j] = 0;
	} else {
	  if (trans_probs[i+nstates*j] > trans_probs_max[j]) {
	    trans_probs_whichmax[j] = i;
	    trans_probs_max[j] = trans_probs[i+nstates*j];
	  }
	}
      }
    }
    vsum = 0;
    for (i=0; i<nstates; i++) {
      v[i+t*nstates] = emit_probs[i+t*nstates] * trans_probs_max[i];
      vsum = vsum + v[i+t*nstates];
    }
    if (vsum > 0) {
      for (i=0; i<nstates; i++) {
	v[i+t*nstates] = v[i+t*nstates]/vsum;
	v_path[i+t*nstates] = trans_probs_whichmax[i];
      }
    } else {
      for (i=0; i<nstates; i++) {
	v[i+t*nstates] = v[i+(t-1)*nstates];
	v_path[i+t*nstates] = v_path[i+(t-1)*nstates];
      }
    }
  }

  for (i=0; i<nstates; i++) {
    if (i==0) {
      path[tmax-1] = 0;
      running_max = v[i+(tmax-1)*nstates];
    } else {
      if (v[i+(tmax-1)*nstates] > running_max) {
	path[tmax-1] = i;
	running_max = v[i+(tmax-1)*nstates];
      }
    }
  }
  for (t=(tmax-2); t>=0; t--) {
    R_CheckUserInterrupt();
    path[t] = v_path[path[t+1]+(t+1)*nstates];
  }
}

static R_NativePrimitiveArgType negloglike_t[] = {INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP};

static R_NativePrimitiveArgType viterbi_t[] = {INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, INTSXP};

static const R_CMethodDef cMethods[] = {
  {"negloglike", (DL_FUNC) &negloglike, 8, negloglike_t},
  {"viterbi", (DL_FUNC) &viterbi, 11, viterbi_t},
  {NULL, NULL, 0}
};

void R_init_exomeCopy(DllInfo *info)
{
  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
}

/*
 * dynamics_generate.c
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "dynamics_generate".
 *
 * Model version              : 1.8
 * Simulink Coder version : 9.7 (R2022a) 13-Nov-2021
 * C source code generated on : Sat Dec 10 17:25:09 2022
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "dynamics_generate.h"
#include <string.h>
#include <stdio.h>
#include "rtwtypes.h"
#include "dynamics_generate_private.h"
#include "rt_nonfinite.h"
#include <stdlib.h>
#include <errno.h>

#define MAX_DIGITS 12

/* Block signals (default storage) */
B_dynamics_generate_T dynamics_generate_B;

/* Continuous states */
X_dynamics_generate_T dynamics_generate_X;

/* Block states (default storage) */
DW_dynamics_generate_T dynamics_generate_DW;

/* External outputs (root outports fed by signals with default storage) */
ExtY_dynamics_generate_T dynamics_generate_Y;

/* Real-time model */
static RT_MODEL_dynamics_generate_T dynamics_generate_M_;
RT_MODEL_dynamics_generate_T *const dynamics_generate_M = &dynamics_generate_M_;

//function definitions
int get_u_input(char* uval);
int get_h_input(char* hval);
void getIntegerFromStdin(int *inputInteger);

/*
 * This function updates continuous states using the ODE3 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  /* Solver Matrices */
  static const real_T rt_ODE3_A[3] = {
    1.0/2.0, 3.0/4.0, 1.0
  };

  static const real_T rt_ODE3_B[3][3] = {
    { 1.0/2.0, 0.0, 0.0 },

    { 0.0, 3.0/4.0, 0.0 },

    { 2.0/9.0, 1.0/3.0, 4.0/9.0 }
  };

  time_T t = rtsiGetT(si);
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE3_IntgData *id = (ODE3_IntgData *)rtsiGetSolverData(si);
  real_T *y = id->y;
  real_T *f0 = id->f[0];
  real_T *f1 = id->f[1];
  real_T *f2 = id->f[2];
  real_T hB[3];
  int_T i;
  int_T nXc = 4;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) memcpy(y, x,
                (uint_T)nXc*sizeof(real_T));

  /* Assumes that rtsiSetT and ModelOutputs are up-to-date */
  /* f0 = f(t,y) */
  rtsiSetdX(si, f0);
  dynamics_generate_derivatives();

  /* f(:,2) = feval(odefile, t + hA(1), y + f*hB(:,1), args(:)(*)); */
  hB[0] = h * rt_ODE3_B[0][0];
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[0]);
  rtsiSetdX(si, f1);
  dynamics_generate_step();
  dynamics_generate_derivatives();

  /* f(:,3) = feval(odefile, t + hA(2), y + f*hB(:,2), args(:)(*)); */
  for (i = 0; i <= 1; i++) {
    hB[i] = h * rt_ODE3_B[1][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[1]);
  rtsiSetdX(si, f2);
  dynamics_generate_step();
  dynamics_generate_derivatives();

  /* tnew = t + hA(3);
     ynew = y + f*hB(:,3); */
  for (i = 0; i <= 2; i++) {
    hB[i] = h * rt_ODE3_B[2][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1] + f2[i]*hB[2]);
  }

  rtsiSetT(si, tnew);
  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model step function */
void dynamics_generate_step(void)
{
  real_T x[4];
  int8_T hnew[4];
  int8_T ii_data[4];
  int8_T newind_data[4];

  //static variable for the counter timestep
  static int count = 0;
  count++;

  if (rtmIsMajorTimeStep(dynamics_generate_M)) {
    /* set solver stop time */
    if (!(dynamics_generate_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&dynamics_generate_M->solverInfo,
                            ((dynamics_generate_M->Timing.clockTickH0 + 1) *
        dynamics_generate_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&dynamics_generate_M->solverInfo,
                            ((dynamics_generate_M->Timing.clockTick0 + 1) *
        dynamics_generate_M->Timing.stepSize0 +
        dynamics_generate_M->Timing.clockTickH0 *
        dynamics_generate_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(dynamics_generate_M)) {
    dynamics_generate_M->Timing.t[0] = rtsiGetT(&dynamics_generate_M->solverInfo);
  }

        //outputs are below (not all are outputs, just the ones that say Outport)
          //print them out here

  /* Integrator: '<S1>/x1 dynamics' */
  dynamics_generate_B.x1 = dynamics_generate_X.x1dynamics_CSTATE;

  /* Outport: '<Root>/x1' */
  dynamics_generate_Y.x1 = dynamics_generate_B.x1;

  /* Integrator: '<S1>/x2 dynamics' */
  dynamics_generate_B.x2 = dynamics_generate_X.x2dynamics_CSTATE;

  /* Outport: '<Root>/x2' */
  dynamics_generate_Y.x2 = dynamics_generate_B.x2;

  /* Integrator: '<S1>/x3 dynamics' */
  dynamics_generate_B.x3 = dynamics_generate_X.x3dynamics_CSTATE;

  /* Outport: '<Root>/x3' */
  dynamics_generate_Y.x3 = dynamics_generate_B.x3;

  /* Integrator: '<S1>/x4 dynamics' */
  dynamics_generate_B.x4 = dynamics_generate_X.x4dynamics_CSTATE;

  /* Outport: '<Root>/x4' */
  dynamics_generate_Y.x4 = dynamics_generate_B.x4;

  //**************
  //print outputs here
  //**************

  printf("**********TIMESTEP %d********\n", count);
  //std::cout << "**********TIMESTEP**********" << std::endl;
  printf("Room 1 Temp: %f\n", dynamics_generate_Y.x1);
  //std::cout << "Room 1 Temp: " << dynamics_generate_Y.x1 << std::endl;
  printf("Room 2 Temp: %f\n", dynamics_generate_Y.x2);
  //std::cout << "Room 2 Temp: " << dynamics_generate_Y.x2 << std::endl;
  printf("Room 3 Temp: %f\n", dynamics_generate_Y.x3);
  //std::cout << "Room 3 Temp: " << dynamics_generate_Y.x3 << std::endl;
  printf("Room 4 Temp: %f\n", dynamics_generate_Y.x4);
  //std::cout << "Room 4 Temp: " << dynamics_generate_Y.x4 << std::endl;
  printf("****************************\n");
  //printf("uval: %f\n", dynamics_generate_P.u);
  //std::cout << "****************************" << std::endl;

  //end print outputs

  if (rtmIsMajorTimeStep(dynamics_generate_M)) {
    real_T biggestdiff;
    real_T curdiff;
    int32_T b_ii;
    int32_T idx;
    boolean_T exitg1;

    /* Memory: '<Root>/Memory' */
    dynamics_generate_B.h1 = dynamics_generate_DW.Memory_PreviousInput;

    /* Memory: '<Root>/Memory1' */
    dynamics_generate_B.h2 = dynamics_generate_DW.Memory1_PreviousInput;

    /* Memory: '<Root>/Memory2' */
    dynamics_generate_B.h3 = dynamics_generate_DW.Memory2_PreviousInput;

    /* Memory: '<Root>/Memory3' */
    dynamics_generate_B.h4 = dynamics_generate_DW.Memory3_PreviousInput;

    /* MATLAB Function: '<S2>/MATLAB Function' incorporates:
     *  Constant: '<Root>/Constant4'
     *  Constant: '<Root>/Constant5'
     *  Constant: '<Root>/Constant6'
     *  Constant: '<Root>/Constant7'
     */
    x[0] = dynamics_generate_B.x1;
    x[1] = dynamics_generate_B.x2;
    x[2] = dynamics_generate_B.x3;
    x[3] = dynamics_generate_B.x4;
    hnew[0] = 0;
    hnew[1] = 0;
    hnew[2] = 0;
    hnew[3] = 0;
    if (!dynamics_generate_DW.locations_not_empty) {      //set h in dynamics_generate_B from input
      dynamics_generate_DW.locations[0] = dynamics_generate_B.h1;
      dynamics_generate_DW.locations[1] = dynamics_generate_B.h2;
      dynamics_generate_DW.locations[2] = dynamics_generate_B.h3;
      dynamics_generate_DW.locations[3] = dynamics_generate_B.h4;
      dynamics_generate_DW.locations_not_empty = true;
    }

    idx = 0;
    b_ii = 0;
    exitg1 = false;
    while ((!exitg1) && (b_ii < 4)) {
      if (dynamics_generate_DW.locations[b_ii] == 1.0) {
        idx++;
        ii_data[idx - 1] = (int8_T)(b_ii + 1);
        if (idx >= 4) {
          exitg1 = true;
        } else {
          b_ii++;
        }
      } else {
        b_ii++;
      }
    }

    if (idx < 1) {
      idx = 0;
    }

    biggestdiff = 0.0;
    idx--;
    if (idx >= 0) {
      memcpy(&newind_data[0], &ii_data[0], (idx + 1) * sizeof(int8_T));
    }

    if ((dynamics_generate_B.x1 <= dynamics_generate_P.get[0]) && (x[ii_data[0]
         - 1] - dynamics_generate_B.x1 >= dynamics_generate_P.dif[0])) {
      curdiff = dynamics_generate_P.on[0] - dynamics_generate_B.x1;
      if (curdiff > 0.0) {
        biggestdiff = curdiff;
        newind_data[0] = 1;
      }
    }

    if ((dynamics_generate_B.x2 <= dynamics_generate_P.get[1]) && (x[ii_data[0]
         - 1] - dynamics_generate_B.x2 >= dynamics_generate_P.dif[1])) {
      curdiff = dynamics_generate_P.on[1] - dynamics_generate_B.x2;
      if (curdiff > biggestdiff) {
        biggestdiff = curdiff;
        newind_data[0] = 2;
      }
    }

    if ((dynamics_generate_B.x3 <= dynamics_generate_P.get[2]) && (x[ii_data[0]
         - 1] - dynamics_generate_B.x3 >= dynamics_generate_P.dif[2])) {
      curdiff = dynamics_generate_P.on[2] - dynamics_generate_B.x3;
      if (curdiff > biggestdiff) {
        biggestdiff = curdiff;
        newind_data[0] = 3;
      }
    }

    if ((dynamics_generate_B.x4 <= dynamics_generate_P.get[3]) && (x[ii_data[0]
         - 1] - dynamics_generate_B.x4 >= dynamics_generate_P.dif[3]) &&
        (dynamics_generate_P.on[3] - dynamics_generate_B.x4 > biggestdiff)) {
      newind_data[0] = 4;
    }

    b_ii = newind_data[0];
    dynamics_generate_DW.locations[b_ii - 1] = 1.0;
    biggestdiff = x[b_ii - 1];
    if ((biggestdiff <= dynamics_generate_P.on[b_ii - 1]) && (biggestdiff <
         dynamics_generate_P.off[b_ii - 1])) {
      hnew[b_ii - 1] = 1;
    }

    biggestdiff = 0.0;
    if ((dynamics_generate_B.x1 <= dynamics_generate_P.get[0]) && (x[ii_data[0]
         - 1] - dynamics_generate_B.x1 >= dynamics_generate_P.dif[0]) && (b_ii
         != 1)) {
      curdiff = dynamics_generate_P.on[0] - dynamics_generate_B.x1;
      if (curdiff > 0.0) {
        biggestdiff = curdiff;
        newind_data[1] = 1;
      }
    }

    if ((dynamics_generate_B.x2 <= dynamics_generate_P.get[1]) && (x[ii_data[0]
         - 1] - dynamics_generate_B.x2 >= dynamics_generate_P.dif[1]) &&
        (newind_data[0] != 2)) {
      curdiff = dynamics_generate_P.on[1] - dynamics_generate_B.x2;
      if (curdiff > biggestdiff) {
        biggestdiff = curdiff;
        newind_data[1] = 2;
      }
    }

    if ((dynamics_generate_B.x3 <= dynamics_generate_P.get[2]) && (x[ii_data[0]
         - 1] - dynamics_generate_B.x3 >= dynamics_generate_P.dif[2]) &&
        (newind_data[0] != 3)) {
      curdiff = dynamics_generate_P.on[2] - dynamics_generate_B.x3;
      if (curdiff > biggestdiff) {
        biggestdiff = curdiff;
        newind_data[1] = 3;
      }
    }

    if ((dynamics_generate_B.x4 <= dynamics_generate_P.get[3]) && (x[ii_data[0]
         - 1] - dynamics_generate_B.x4 >= dynamics_generate_P.dif[3]) &&
        (newind_data[0] != 4) && (dynamics_generate_P.on[3] -
         dynamics_generate_B.x4 > biggestdiff)) {
      newind_data[1] = 4;
    }

    idx = newind_data[1] - 1;
    dynamics_generate_DW.locations[idx] = 1.0;
    biggestdiff = x[idx];
    if ((biggestdiff <= dynamics_generate_P.on[idx]) && (biggestdiff <
         dynamics_generate_P.off[idx])) {
      hnew[idx] = 1;
    }

    dynamics_generate_B.h1 = hnew[0];
    dynamics_generate_B.h2 = hnew[1];
    dynamics_generate_B.h3 = hnew[2];
    dynamics_generate_B.h4 = hnew[3];

    /* End of MATLAB Function: '<S2>/MATLAB Function' */

    /* Gain: '<S1>/Gain' */
    dynamics_generate_B.Gain = dynamics_generate_P.c[0] * dynamics_generate_B.h1;

    /* Gain: '<S1>/c2' */
    dynamics_generate_B.c2 = dynamics_generate_P.c[1] * dynamics_generate_B.h2;

    /* Gain: '<S1>/c3' */
    dynamics_generate_B.c3 = dynamics_generate_P.c[2] * dynamics_generate_B.h3;

    /* Gain: '<S1>/c4' */
    dynamics_generate_B.c4 = dynamics_generate_P.c[3] * dynamics_generate_B.h4;
  }

  /* Sum: '<S1>/Add' incorporates:
   *  Constant: '<S1>/u'
   *  Gain: '<S1>/a21'
   *  Gain: '<S1>/a23'
   *  Gain: '<S1>/a24'
   *  Gain: '<S1>/b2'
   *  Sum: '<S1>/(2) b + c'
   *  Sum: '<S1>/2 Summation'
   *  Sum: '<S1>/u-x2'
   *  Sum: '<S1>/x1-x2'
   *  Sum: '<S1>/x3-x2'
   *  Sum: '<S1>/x4-x2'
   */
  dynamics_generate_B.x2_dot = (((dynamics_generate_B.x1 -
    dynamics_generate_B.x2) * dynamics_generate_P.a21_Gain +
    (dynamics_generate_B.x3 - dynamics_generate_B.x2) *
    dynamics_generate_P.a23_Gain) + (dynamics_generate_B.x4 -
    dynamics_generate_B.x2) * dynamics_generate_P.a24_Gain) +
    ((dynamics_generate_P.u - dynamics_generate_B.x2) * dynamics_generate_P.b[1]
     + dynamics_generate_B.c2);

  /* Sum: '<S1>/Add1' incorporates:
   *  Constant: '<S1>/u'
   *  Gain: '<S1>/a12'
   *  Gain: '<S1>/a13'
   *  Gain: '<S1>/a14'
   *  Gain: '<S1>/b1'
   *  Sum: '<S1>/(1) b + c'
   *  Sum: '<S1>/1 Summation'
   *  Sum: '<S1>/u-x1'
   *  Sum: '<S1>/x2-x1'
   *  Sum: '<S1>/x3-x1'
   *  Sum: '<S1>/x4-x1'
   */
  dynamics_generate_B.x1_dot = (((dynamics_generate_B.x2 -
    dynamics_generate_B.x1) * dynamics_generate_P.a12_Gain +
    (dynamics_generate_B.x3 - dynamics_generate_B.x1) *
    dynamics_generate_P.a13_Gain) + (dynamics_generate_B.x4 -
    dynamics_generate_B.x1) * dynamics_generate_P.a14_Gain) +
    ((dynamics_generate_P.u - dynamics_generate_B.x1) * dynamics_generate_P.b[0]
     + dynamics_generate_B.Gain);

  /* Sum: '<S1>/Add2' incorporates:
   *  Constant: '<S1>/u'
   *  Gain: '<S1>/a31'
   *  Gain: '<S1>/a32'
   *  Gain: '<S1>/a34'
   *  Gain: '<S1>/b3'
   *  Sum: '<S1>/(3) b + c'
   *  Sum: '<S1>/3 Summation'
   *  Sum: '<S1>/u - x3'
   *  Sum: '<S1>/x1-x3'
   *  Sum: '<S1>/x2-x3'
   *  Sum: '<S1>/x4-x3'
   */
  dynamics_generate_B.x3_dot = (((dynamics_generate_B.x1 -
    dynamics_generate_B.x3) * dynamics_generate_P.a31_Gain +
    (dynamics_generate_B.x2 - dynamics_generate_B.x3) *
    dynamics_generate_P.a32_Gain) + (dynamics_generate_B.x4 -
    dynamics_generate_B.x3) * dynamics_generate_P.a34_Gain) +
    ((dynamics_generate_P.u - dynamics_generate_B.x3) * dynamics_generate_P.b[2]
     + dynamics_generate_B.c3);

  /* Sum: '<S1>/Add3' incorporates:
   *  Constant: '<S1>/u'
   *  Gain: '<S1>/Gain1'
   *  Gain: '<S1>/a41'
   *  Gain: '<S1>/a42'
   *  Gain: '<S1>/a43'
   *  Sum: '<S1>/(4) b + c'
   *  Sum: '<S1>/4 Summation term'
   *  Sum: '<S1>/u - x4'
   *  Sum: '<S1>/x1-x4'
   *  Sum: '<S1>/x2-x4'
   *  Sum: '<S1>/x3-x4'
   */
  dynamics_generate_B.x4_dot = (((dynamics_generate_B.x1 -
    dynamics_generate_B.x4) * dynamics_generate_P.a41_Gain +
    (dynamics_generate_B.x2 - dynamics_generate_B.x4) *
    dynamics_generate_P.a42_Gain) + (dynamics_generate_B.x3 -
    dynamics_generate_B.x4) * dynamics_generate_P.a43_Gain) +
    ((dynamics_generate_P.u - dynamics_generate_B.x4) * dynamics_generate_P.b[3]
     + dynamics_generate_B.c4);
  if (rtmIsMajorTimeStep(dynamics_generate_M)) {
    /* Matfile logging */
    rt_UpdateTXYLogVars(dynamics_generate_M->rtwLogInfo,
                        (dynamics_generate_M->Timing.t));
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(dynamics_generate_M)) {
    if (rtmIsMajorTimeStep(dynamics_generate_M)) {
      /* Update for Memory: '<Root>/Memory' */
      dynamics_generate_DW.Memory_PreviousInput = dynamics_generate_B.h1;

      /* Update for Memory: '<Root>/Memory1' */
      dynamics_generate_DW.Memory1_PreviousInput = dynamics_generate_B.h2;

      /* Update for Memory: '<Root>/Memory2' */
      dynamics_generate_DW.Memory2_PreviousInput = dynamics_generate_B.h3;

      /* Update for Memory: '<Root>/Memory3' */
      dynamics_generate_DW.Memory3_PreviousInput = dynamics_generate_B.h4;
    }
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(dynamics_generate_M)) {
    /* signal main to stop simulation */
    {                                  /* Sample time: [0.0s, 0.0s] */
      if ((rtmGetTFinal(dynamics_generate_M)!=-1) &&
          !((rtmGetTFinal(dynamics_generate_M)-
             (((dynamics_generate_M->Timing.clockTick1+
                dynamics_generate_M->Timing.clockTickH1* 4294967296.0)) * 0.1)) >
            (((dynamics_generate_M->Timing.clockTick1+
               dynamics_generate_M->Timing.clockTickH1* 4294967296.0)) * 0.1) *
            (DBL_EPSILON))) {
        rtmSetErrorStatus(dynamics_generate_M, "Simulation finished");
      }
    }

    rt_ertODEUpdateContinuousStates(&dynamics_generate_M->solverInfo);

    /* Update absolute time for base rate */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick0 and the high bits
     * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++dynamics_generate_M->Timing.clockTick0)) {
      ++dynamics_generate_M->Timing.clockTickH0;
    }

    dynamics_generate_M->Timing.t[0] = rtsiGetSolverStopTime
      (&dynamics_generate_M->solverInfo);

    {
      /* Update absolute timer for sample time: [0.1s, 0.0s] */
      /* The "clockTick1" counts the number of times the code of this task has
       * been executed. The resolution of this integer timer is 0.1, which is the step size
       * of the task. Size of "clockTick1" ensures timer will not overflow during the
       * application lifespan selected.
       * Timer of this task consists of two 32 bit unsigned integers.
       * The two integers represent the low bits Timing.clockTick1 and the high bits
       * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
       */
      dynamics_generate_M->Timing.clockTick1++;
      if (!dynamics_generate_M->Timing.clockTick1) {
        dynamics_generate_M->Timing.clockTickH1++;
      }
    }
  }                                    /* end MajorTimeStep */
}

/* Derivatives for root system: '<Root>' */
void dynamics_generate_derivatives(void)
{
  XDot_dynamics_generate_T *_rtXdot;
  _rtXdot = ((XDot_dynamics_generate_T *) dynamics_generate_M->derivs);

  /* Derivatives for Integrator: '<S1>/x1 dynamics' */
  _rtXdot->x1dynamics_CSTATE = dynamics_generate_B.x1_dot;

  /* Derivatives for Integrator: '<S1>/x2 dynamics' */
  _rtXdot->x2dynamics_CSTATE = dynamics_generate_B.x2_dot;

  /* Derivatives for Integrator: '<S1>/x3 dynamics' */
  _rtXdot->x3dynamics_CSTATE = dynamics_generate_B.x3_dot;

  /* Derivatives for Integrator: '<S1>/x4 dynamics' */
  _rtXdot->x4dynamics_CSTATE = dynamics_generate_B.x4_dot;
}

/* Model initialize function */
void dynamics_generate_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* initialize real-time model */
  (void) memset((void *)dynamics_generate_M, 0,
                sizeof(RT_MODEL_dynamics_generate_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&dynamics_generate_M->solverInfo,
                          &dynamics_generate_M->Timing.simTimeStep);
    rtsiSetTPtr(&dynamics_generate_M->solverInfo, &rtmGetTPtr
                (dynamics_generate_M));
    rtsiSetStepSizePtr(&dynamics_generate_M->solverInfo,
                       &dynamics_generate_M->Timing.stepSize0);
    rtsiSetdXPtr(&dynamics_generate_M->solverInfo, &dynamics_generate_M->derivs);
    rtsiSetContStatesPtr(&dynamics_generate_M->solverInfo, (real_T **)
                         &dynamics_generate_M->contStates);
    rtsiSetNumContStatesPtr(&dynamics_generate_M->solverInfo,
      &dynamics_generate_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&dynamics_generate_M->solverInfo,
      &dynamics_generate_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&dynamics_generate_M->solverInfo,
      &dynamics_generate_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&dynamics_generate_M->solverInfo,
      &dynamics_generate_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&dynamics_generate_M->solverInfo, (&rtmGetErrorStatus
      (dynamics_generate_M)));
    rtsiSetRTModelPtr(&dynamics_generate_M->solverInfo, dynamics_generate_M);
  }

  rtsiSetSimTimeStep(&dynamics_generate_M->solverInfo, MAJOR_TIME_STEP);
  dynamics_generate_M->intgData.y = dynamics_generate_M->odeY;
  dynamics_generate_M->intgData.f[0] = dynamics_generate_M->odeF[0];
  dynamics_generate_M->intgData.f[1] = dynamics_generate_M->odeF[1];
  dynamics_generate_M->intgData.f[2] = dynamics_generate_M->odeF[2];
  dynamics_generate_M->contStates = ((X_dynamics_generate_T *)
    &dynamics_generate_X);
  rtsiSetSolverData(&dynamics_generate_M->solverInfo, (void *)
                    &dynamics_generate_M->intgData);
  rtsiSetIsMinorTimeStepWithModeChange(&dynamics_generate_M->solverInfo, false);
  rtsiSetSolverName(&dynamics_generate_M->solverInfo,"ode3");
  rtmSetTPtr(dynamics_generate_M, &dynamics_generate_M->Timing.tArray[0]);
  rtmSetTFinal(dynamics_generate_M, 100.0);
  dynamics_generate_M->Timing.stepSize0 = 0.1;

  /* Setup for data logging */
  {
    static RTWLogInfo rt_DataLoggingInfo;
    rt_DataLoggingInfo.loggingInterval = (NULL);
    dynamics_generate_M->rtwLogInfo = &rt_DataLoggingInfo;
  }

  /* Setup for data logging */
  {
    rtliSetLogXSignalInfo(dynamics_generate_M->rtwLogInfo, (NULL));
    rtliSetLogXSignalPtrs(dynamics_generate_M->rtwLogInfo, (NULL));
    rtliSetLogT(dynamics_generate_M->rtwLogInfo, "tout");
    rtliSetLogX(dynamics_generate_M->rtwLogInfo, "");
    rtliSetLogXFinal(dynamics_generate_M->rtwLogInfo, "");
    rtliSetLogVarNameModifier(dynamics_generate_M->rtwLogInfo, "rt_");
    rtliSetLogFormat(dynamics_generate_M->rtwLogInfo, 4);
    rtliSetLogMaxRows(dynamics_generate_M->rtwLogInfo, 0);
    rtliSetLogDecimation(dynamics_generate_M->rtwLogInfo, 1);
    rtliSetLogY(dynamics_generate_M->rtwLogInfo, "");
    rtliSetLogYSignalInfo(dynamics_generate_M->rtwLogInfo, (NULL));
    rtliSetLogYSignalPtrs(dynamics_generate_M->rtwLogInfo, (NULL));
  }

  /* block I/O */
  (void) memset(((void *) &dynamics_generate_B), 0,
                sizeof(B_dynamics_generate_T));

  /* states (continuous) */
  {
    (void) memset((void *)&dynamics_generate_X, 0,
                  sizeof(X_dynamics_generate_T));
  }

  /* states (dwork) */
  (void) memset((void *)&dynamics_generate_DW, 0,
                sizeof(DW_dynamics_generate_T));

  /* external outputs */
  (void)memset(&dynamics_generate_Y, 0, sizeof(ExtY_dynamics_generate_T));

  /* Matfile logging */
  rt_StartDataLoggingWithStartTime(dynamics_generate_M->rtwLogInfo, 0.0,
    rtmGetTFinal(dynamics_generate_M), dynamics_generate_M->Timing.stepSize0,
    (&rtmGetErrorStatus(dynamics_generate_M)));

  /* InitializeConditions for Integrator: '<S1>/x1 dynamics' */
  dynamics_generate_X.x1dynamics_CSTATE = dynamics_generate_P.x0[0];

  /* InitializeConditions for Integrator: '<S1>/x2 dynamics' */
  dynamics_generate_X.x2dynamics_CSTATE = dynamics_generate_P.x0[1];

  /* InitializeConditions for Integrator: '<S1>/x3 dynamics' */
  dynamics_generate_X.x3dynamics_CSTATE = dynamics_generate_P.x0[2];

  /* InitializeConditions for Integrator: '<S1>/x4 dynamics' */
  dynamics_generate_X.x4dynamics_CSTATE = dynamics_generate_P.x0[3];

  //**********
  //looks like this is where we set h for the first time
    //for u, look for dynamics_generate_P.u --> likely just replace dynamics_generate_P.u
    //with input from the file as well
  //now, just need to output to std::cout when we calculate x1, x2, x3, x4
    //what should the output look like? 
    //TIME: 1
      //Room1 Temp: 15.3
      //Room2 Temp: 16.2
      //...
    //TIME: 2
  //***********

  //collect inputs here
  char* h1str = "h1";
  char* h2str = "h2";
  char* h3str = "h3";
  char* h4str = "h4";
  char* ustr = "u";
  printf("***Enter in the necessary inputs for the executable.\n");
  printf("***Inputs are h1, h2, h3, h4, and u\n");
  printf("***h1, h2, h3, and h4 are all integers either 0 or 1. Their sum must be exactly 2.");
  printf("***u is an integer with any value. Note that the problem was defined initially with u=6.\n");
  printf("***If u input deviates too much from this value, especially if it is too large, the system may not converge.\n");
  printf("***(Initial definitions: u = 6, h1 = 0, h2 = 1, h3 = 1, h4 = 0.)\n");
  int h1 = get_h_input(h1str);
  printf("h1 value entered: %d\n", h1);
  int h2 = get_h_input(h2str);
  printf("h2 value entered: %d\n", h2);
  int h3 = get_h_input(h1str);
  printf("h3 value entered: %d\n", h3);
  int h4 = get_h_input(h2str);
  printf("h4 value entered: %d\n", h4);
  int u = get_u_input(ustr);
  printf("u value entered: %d\n", u);
  dynamics_generate_P.u = u;

  /* InitializeConditions for Memory: '<Root>/Memory' */
  dynamics_generate_DW.Memory_PreviousInput = dynamics_generate_P.h[0];

  /* InitializeConditions for Memory: '<Root>/Memory1' */
  dynamics_generate_DW.Memory1_PreviousInput = dynamics_generate_P.h[1];

  /* InitializeConditions for Memory: '<Root>/Memory2' */
  dynamics_generate_DW.Memory2_PreviousInput = dynamics_generate_P.h[2];

  /* InitializeConditions for Memory: '<Root>/Memory3' */
  dynamics_generate_DW.Memory3_PreviousInput = dynamics_generate_P.h[3];

  /* SystemInitialize for MATLAB Function: '<S2>/MATLAB Function' */
  dynamics_generate_DW.locations_not_empty = false;
}

/* Model terminate function */
void dynamics_generate_terminate(void)
{
  /* (no terminate code required) */
}

int get_h_input(char* hval) {
  int h = 2;
  while (h != 1 && h != 0) {
    printf("Enter a value for %s. Make sure to enter either 0 or 1: ", hval);
    getIntegerFromStdin(&h);
  }
  return h;
}

int get_u_input(char* uval) {
  int u;
  printf("Enter a value for %s. Make sure to enter an integer: ", uval);
  getIntegerFromStdin(&u);
  return u;
}

//the following code was taken from:
  //https://dev-notes.eu/2019/05/Integer-Input-in-C/#:~:text=To%20scan%20integer%20input%2C%20first,not%20be%20the%20best%20choice.
static inline void ClearInputBuffer() 
{
	char c = 0;
	// Loop over input buffer and consume chars until buffer is empty
	while ((c = getchar()) != '\n' && c != EOF);
}

void getIntegerFromStdin(int *inputInteger)
{
	char *inputBuffer = malloc(sizeof(char) * MAX_DIGITS);
	memset(inputBuffer, 0, MAX_DIGITS);
	char *input = NULL;
	while (input == NULL) {
		// Note that fgets returns inputBuffer on success.
		// This becomes important when freeing - free either `input` or
		// `inputBuffer` to avoid an attempted double-free error.
		input = fgets(inputBuffer, MAX_DIGITS, stdin);
		
		// If fgets() receives less than MAX_DIGITS, the last char in the array is '\n'.
		// Therefore if the last char is not '\n', too many characters were entered.
		if (inputBuffer[strlen(inputBuffer) - 1] != '\n') {
			fprintf(stderr, "[ERROR]: Too many characters: max input is %d chars.\n", MAX_DIGITS);
			ClearInputBuffer();
			input = NULL;
			continue;
		}

		// Check that the input can be intepreted as an integer
		// Convert to integer using `strtol()`
		errno = 0;
		char *endptr = NULL;
		*inputInteger = strtol(input, &endptr, 10);
		
		// If an integer was not found, endptr remains set to input
		if (input == endptr) {
			// Remove trailing newline by adding NUL at the index of the
			// terminating '\n' character. See man strcspn - this function
			// gets the length of a prefix substring.
			input[strcspn(input, "\n")] = 0;
			printf("Invalid input: no integer found in %s.\n", input);
			input = NULL;
		}
		if (errno != 0) {
			fprintf(stderr, "[ERROR]: That doesn't look like an integer.\n");
			input = NULL;
		}
	}
	free(inputBuffer);
}
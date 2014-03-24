/*
 * MATLAB Compiler: 3.0
 * Date: Thu Aug 28 15:51:06 2003
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-B" "sgl" "-m" "-W"
 * "main" "-L" "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "-W" "mainhg"
 * "libmwsglm.mlib" "process_rowdata" 
 */
#include "detectstepsplot.h"
#include "khazoomsafe.h"
#include "libmatlbm.h"
#include "libmmfile.h"
#include "libmwsglm.h"
#include "movelinex.h"

static mxChar _array1_[2] = { 'o', 'n' };
static mxArray * _mxarray0_;
static mxArray * _mxarray2_;

static mxChar _array4_[1] = { 'g' };
static mxArray * _mxarray3_;

static mxChar _array6_[9] = { 'E', 'r', 'a', 's', 'e', 'M', 'o', 'd', 'e' };
static mxArray * _mxarray5_;

static mxChar _array8_[3] = { 'x', 'o', 'r' };
static mxArray * _mxarray7_;

static mxChar _array10_[5] = { 'Z', 'D', 'a', 't', 'a' };
static mxArray * _mxarray9_;

static double _array12_[2] = { 1.0, 1.0 };
static mxArray * _mxarray11_;

static mxChar _array16_[67] = { 'D', 'r', 'a', 'g', ' ', 't', 'h', 'e', ' ',
                                'g', 'r', 'e', 'e', 'n', ' ', 'l', 'i', 'n',
                                'e', 's', ' ', 'i', 'n', ' ', 't', 'h', 'e',
                                ' ', 'p', 'l', 'o', 't', ' ', 't', 'o', ' ',
                                'a', 'd', 'j', 'u', 's', 't', ' ', 't', 'h',
                                'e', ' ', 's', 't', 'a', 'r', 't', ' ', 'o',
                                'f', ' ', 't', 'h', 'e', ' ', 'c', 'y', 'c',
                                'l', 'e', 's', '.' };
static mxArray * _mxarray15_;

static mxChar _array18_[95] = { 'A', 'd', 'd', ' ', 'l', 'i', 'n', 'e', 's',
                                ' ', 'b', 'y', ' ', 'r', 'i', 'g', 'h', 't',
                                ' ', 'c', 'l', 'i', 'c', 'k', 'i', 'n', 'g',
                                ' ', 'i', 'n', ' ', 't', 'h', 'e', ' ', 'f',
                                'i', 'g', 'u', 'r', 'e', '.', ' ', 'R', 'e',
                                'm', 'o', 'v', 'e', ' ', 'l', 'i', 'n', 'e',
                                's', ' ', 'b', 'y', ' ', 'h', 'o', 'l', 'd',
                                'i', 'n', 'g', ' ', 'd', 'o', 'w', 'n', ' ',
                                0x0027, 's', 'h', 'i', 'f', 't', 0x0027,
                                ' ', 'w', 'h', 'i', 'l', 'e', ' ', 'c', 'l',
                                'i', 'c', 'k', 'i', 'n', 'g', '.' };
static mxArray * _mxarray17_;

static mxChar _array20_[84] = { 'Z', 'o', 'o', 'm', ' ', 'b', 'y', ' ', 'd',
                                'r', 'a', 'w', 'i', 'n', 'g', ' ', 'o', 'v',
                                'e', 'r', ' ', 't', 'h', 'e', ' ', 'r', 'e',
                                'g', 'i', 'o', 'n', ' ', 'o', 'f', ' ', ' ',
                                'i', 'n', 't', 'e', 'r', 'e', 's', 't', '.',
                                ' ', 'D', 'o', 'u', 'b', 'l', 'e', '-', 'c',
                                'l', 'i', 'c', 'k', ' ', 't', 'o', ' ', 'z',
                                'o', 'o', 'm', ' ', 'b', 'a', 'c', 'k', ' ',
                                't', 'o', ' ', 'o', 'r', 'i', 'g', 'i', 'n',
                                'a', 'l', '.' };
static mxArray * _mxarray19_;

static mxChar _array22_[22] = { 'C', 'l', 'i', 'c', 'k', ' ', 'o', 'k',
                                ' ', 'w', 'h', 'e', 'n', ' ', 'f', 'i',
                                'n', 'i', 's', 'h', 'e', 'd' };
static mxArray * _mxarray21_;

static mxArray * _array14_[4] = { NULL /*_mxarray15_*/, NULL /*_mxarray17_*/,
                                  NULL /*_mxarray19_*/, NULL /*_mxarray21_*/ };
static mxArray * _mxarray13_;

static mxChar _array26_[65] = { 'D', 'r', 'a', 'g', ' ', 't', 'h', 'e', ' ',
                                'g', 'r', 'e', 'e', 'n', ' ', 'l', 'i', 'n',
                                'e', ' ', 'i', 'n', ' ', 't', 'h', 'e', ' ',
                                'p', 'l', 'o', 't', ' ', 't', 'o', ' ', 'a',
                                'd', 'j', 'u', 's', 't', ' ', 't', 'h', 'e',
                                ' ', 's', 't', 'a', 'r', 't', ' ', 'o', 'f',
                                ' ', 't', 'h', 'e', ' ', 'c', 'y', 'c', 'l',
                                'e', '.' };
static mxArray * _mxarray25_;

static mxArray * _array24_[3] = { NULL /*_mxarray25_*/, NULL /*_mxarray19_*/,
                                  NULL /*_mxarray21_*/ };
static mxArray * _mxarray23_;

static mxArray * _array28_[3] = { NULL /*_mxarray17_*/, NULL /*_mxarray19_*/,
                                  NULL /*_mxarray21_*/ };
static mxArray * _mxarray27_;
static mxArray * _mxarray29_;

static mxChar _array31_[12] = { 'I', 'n', 's', 't', 'r', 'u',
                                'c', 't', 'i', 'o', 'n', 's' };
static mxArray * _mxarray30_;

static mxChar _array33_[8] = { 'C', 'h', 'i', 'l', 'd', 'r', 'e', 'n' };
static mxArray * _mxarray32_;
static mxArray * _mxarray34_;

static mxChar _array36_[5] = { 'X', 'D', 'a', 't', 'a' };
static mxArray * _mxarray35_;

void InitializeModule_detectstepsplot(void) {
    _mxarray0_ = mclInitializeString(2, _array1_);
    _mxarray2_ = mclInitializeDouble(1.0);
    _mxarray3_ = mclInitializeString(1, _array4_);
    _mxarray5_ = mclInitializeString(9, _array6_);
    _mxarray7_ = mclInitializeString(3, _array8_);
    _mxarray9_ = mclInitializeString(5, _array10_);
    _mxarray11_ = mclInitializeDoubleVector(2, 1, _array12_);
    _mxarray15_ = mclInitializeString(67, _array16_);
    _array14_[0] = _mxarray15_;
    _mxarray17_ = mclInitializeString(95, _array18_);
    _array14_[1] = _mxarray17_;
    _mxarray19_ = mclInitializeString(84, _array20_);
    _array14_[2] = _mxarray19_;
    _mxarray21_ = mclInitializeString(22, _array22_);
    _array14_[3] = _mxarray21_;
    _mxarray13_ = mclInitializeCellVector(1, 4, _array14_);
    _mxarray25_ = mclInitializeString(65, _array26_);
    _array24_[0] = _mxarray25_;
    _array24_[1] = _mxarray19_;
    _array24_[2] = _mxarray21_;
    _mxarray23_ = mclInitializeCellVector(1, 3, _array24_);
    _array28_[0] = _mxarray17_;
    _array28_[1] = _mxarray19_;
    _array28_[2] = _mxarray21_;
    _mxarray27_ = mclInitializeCellVector(1, 3, _array28_);
    _mxarray29_ = mclInitializeDouble(40.0);
    _mxarray30_ = mclInitializeString(12, _array31_);
    _mxarray32_ = mclInitializeString(8, _array33_);
    _mxarray34_ = mclInitializeDoubleVector(0, 0, (double *)NULL);
    _mxarray35_ = mclInitializeString(5, _array36_);
}

void TerminateModule_detectstepsplot(void) {
    mxDestroyArray(_mxarray35_);
    mxDestroyArray(_mxarray34_);
    mxDestroyArray(_mxarray32_);
    mxDestroyArray(_mxarray30_);
    mxDestroyArray(_mxarray29_);
    mxDestroyArray(_mxarray27_);
    mxDestroyArray(_mxarray23_);
    mxDestroyArray(_mxarray25_);
    mxDestroyArray(_mxarray13_);
    mxDestroyArray(_mxarray21_);
    mxDestroyArray(_mxarray19_);
    mxDestroyArray(_mxarray17_);
    mxDestroyArray(_mxarray15_);
    mxDestroyArray(_mxarray11_);
    mxDestroyArray(_mxarray9_);
    mxDestroyArray(_mxarray7_);
    mxDestroyArray(_mxarray5_);
    mxDestroyArray(_mxarray3_);
    mxDestroyArray(_mxarray2_);
    mxDestroyArray(_mxarray0_);
}

static mxArray * Mdetectstepsplot(mxArray * * figh,
                                  int nargout_,
                                  mxArray * stpfr,
                                  mxArray * mdf);

_mexLocalFunctionTable _local_function_table_detectstepsplot
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfDetectstepsplot" contains the normal interface for the
 * "detectstepsplot" M-function from file
 * "k:\matlab\rowpower\detectstepsplot.m" (lines 1-78). This function processes
 * any input arguments and passes them to the implementation version of the
 * function, appearing above.
 */
mxArray * mlfDetectstepsplot(mxArray * * figh, mxArray * stpfr, mxArray * mdf) {
    int nargout = 1;
    mxArray * cyclefr = NULL;
    mxArray * figh__ = NULL;
    mlfEnterNewContext(1, 2, figh, stpfr, mdf);
    if (figh != NULL) {
        ++nargout;
    }
    cyclefr = Mdetectstepsplot(&figh__, nargout, stpfr, mdf);
    mlfRestorePreviousContext(1, 2, figh, stpfr, mdf);
    if (figh != NULL) {
        mclCopyOutputArg(figh, figh__);
    } else {
        mxDestroyArray(figh__);
    }
    return mlfReturnValue(cyclefr);
}

/*
 * The function "mlxDetectstepsplot" contains the feval interface for the
 * "detectstepsplot" M-function from file
 * "k:\matlab\rowpower\detectstepsplot.m" (lines 1-78). The feval function
 * calls the implementation version of detectstepsplot through this function.
 * This function processes any input arguments and passes them to the
 * implementation version of the function, appearing above.
 */
void mlxDetectstepsplot(int nlhs,
                        mxArray * plhs[],
                        int nrhs,
                        mxArray * prhs[]) {
    mxArray * mprhs[2];
    mxArray * mplhs[2];
    int i;
    if (nlhs > 2) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: detectstepsplot Line: 1 Colum"
            "n: 1 The function \"detectstepsplot\" was called wi"
            "th more than the declared number of outputs (2)."),
          NULL);
    }
    if (nrhs > 2) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: detectstepsplot Line: 1 Colum"
            "n: 1 The function \"detectstepsplot\" was called wi"
            "th more than the declared number of inputs (2)."),
          NULL);
    }
    for (i = 0; i < 2; ++i) {
        mplhs[i] = NULL;
    }
    for (i = 0; i < 2 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 2; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 2, mprhs[0], mprhs[1]);
    mplhs[0] = Mdetectstepsplot(&mplhs[1], nlhs, mprhs[0], mprhs[1]);
    mlfRestorePreviousContext(0, 2, mprhs[0], mprhs[1]);
    plhs[0] = mplhs[0];
    for (i = 1; i < 2 && i < nlhs; ++i) {
        plhs[i] = mplhs[i];
    }
    for (; i < 2; ++i) {
        mxDestroyArray(mplhs[i]);
    }
}

/*
 * The function "Mdetectstepsplot" is the implementation version of the
 * "detectstepsplot" M-function from file
 * "k:\matlab\rowpower\detectstepsplot.m" (lines 1-78). It contains the actual
 * compiled code for that M-function. It is a static function and must only be
 * called from one of the interface functions, appearing below.
 */
/*
 * function  [cyclefr, figh]=detectstepsplot(stpfr,mdf);
 */
static mxArray * Mdetectstepsplot(mxArray * * figh,
                                  int nargout_,
                                  mxArray * stpfr,
                                  mxArray * mdf) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(
          &_local_function_table_detectstepsplot);
    mxArray * cyclefr = NULL;
    mxArray * xd = NULL;
    mxArray * i = NULL;
    mxArray * lines = NULL;
    mxArray * msg = NULL;
    mxArray * cyclelines = NULL;
    mxArray * ans = NULL;
    mxArray * fig = NULL;
    mclCopyArray(&stpfr);
    mclCopyArray(&mdf);
    /*
     * %
     * %  [cyclefr, figh]=detectstepsplot(stepfr,mdf);
     * % 
     * % Plots the time series mdf and over lies the beginning of the cycles as 
     * % vertical lines. The user is then asked to modify the set of lines so 
     * % as to truly correspond to the beginning (and end) of each cycle.
     * 
     * % Kjartan Halvorsen
     * % 2001-02-22
     * 
     * 
     * fig=figure;
     */
    mlfAssign(&fig, mlfNFigure(1, NULL));
    /*
     * figh=fig;
     */
    mlfAssign(figh, mclVv(fig, "fig"));
    /*
     * 
     * plot(mdf);
     */
    mclAssignAns(&ans, mlfNPlot(0, mclVa(mdf, "mdf"), NULL));
    /*
     * hold on
     */
    mlfHold(_mxarray0_);
    /*
     * if (~isempty(stpfr))
     */
    if (mclNotBool(mlfIsempty(mclVa(stpfr, "stpfr")))) {
        /*
         * cyclelines=plot([stpfr';stpfr'],...
         */
        mlfAssign(
          &cyclelines,
          mlfNPlot(
            1,
            mlfVertcat(
              mlfCtranspose(mclVa(stpfr, "stpfr")),
              mlfCtranspose(mclVa(stpfr, "stpfr")),
              NULL),
            mlfVertcat(
              mclMtimes(
                mlfMin(
                  NULL,
                  mlfMin(NULL, mclVa(mdf, "mdf"), NULL, NULL),
                  NULL,
                  NULL),
                mlfOnes(
                  _mxarray2_,
                  mlfScalar(mclLengthInt(mclVa(stpfr, "stpfr"))),
                  NULL)),
              mclMtimes(
                mlfMax(
                  NULL,
                  mlfMax(NULL, mclVa(mdf, "mdf"), NULL, NULL),
                  NULL,
                  NULL),
                mlfOnes(
                  _mxarray2_,
                  mlfScalar(mclLengthInt(mclVa(stpfr, "stpfr"))),
                  NULL)),
              NULL),
            _mxarray3_,
            _mxarray5_,
            _mxarray7_,
            NULL));
        /*
         * [min(min(mdf))*ones(1,length(stpfr))
         * max(max(mdf))*ones(1,length(stpfr))],...
         * 'g','EraseMode','xor');
         * set(cyclelines,'ZData',[1;1]);
         */
        mclAssignAns(
          &ans,
          mlfNSet(
            0, mclVv(cyclelines, "cyclelines"), _mxarray9_, _mxarray11_, NULL));
    /*
     * end
     */
    }
    /*
     * 
     * movelinex(fig);
     */
    mlfMovelinex(mclVv(fig, "fig"), NULL);
    /*
     * 
     * khazoomsafe('on');
     */
    mlfKhazoomsafe(_mxarray0_, NULL);
    /*
     * 
     * if (~isempty(stpfr))
     */
    if (mclNotBool(mlfIsempty(mclVa(stpfr, "stpfr")))) {
        /*
         * if (length(stpfr)>1)
         */
        if (mclLengthInt(mclVa(stpfr, "stpfr")) > 1) {
            /*
             * msg={['Drag the green lines in the plot to adjust ',...
             */
            mlfAssign(&msg, _mxarray13_);
        /*
         * 'the start of the cycles.'],...
         * ['Add lines by right clicking', ...
         * ' in the figure. Remove lines by holding down ''shift''', ...
         * ' while clicking.'],...
         * ['Zoom by drawing over the region of ',...
         * ' interest. Double-click to zoom back to original.'],...
         * 'Click ok when finished'};
         * else
         */
        } else {
            /*
             * msg={['Drag the green line in the plot to adjust ',...
             */
            mlfAssign(&msg, _mxarray23_);
        /*
         * 'the start of the cycle.'],...
         * ['Zoom by drawing over the region of ',...
         * ' interest. Double-click to zoom back to original.'],...
         * 'Click ok when finished'};
         * 
         * end	  
         */
        }
    /*
     * else
     */
    } else {
        /*
         * msg={['Add lines by right clicking', ...
         */
        mlfAssign(&msg, _mxarray27_);
    /*
     * ' in the figure. Remove lines by holding down ''shift''', ...
     * ' while clicking.'],['Zoom by drawing over the region of ',...
     * ' interest. Double-click to zoom back to original.'],...
     * 'Click ok when finished'};
     * end
     */
    }
    /*
     * 
     * msg=textwrap(msg,40);
     */
    mlfAssign(
      &msg,
      mlfNTextwrap(
        1, mlfVarargout(NULL), mclVv(msg, "msg"), _mxarray29_, NULL));
    /*
     * uiwait(msgbox(msg,'Instructions'));
     */
    mlfUiwait(
      mlfNMsgbox(0, mclValueVarargout(), mclVv(msg, "msg"), _mxarray30_, NULL));
    /*
     * 
     * % Check the XData values of the cyclelines
     * figure(fig)
     */
    mclPrintAns(&ans, mlfNFigure(0, mclVv(fig, "fig"), NULL));
    /*
     * lines=get(gca,'Children');
     */
    mlfAssign(&lines, mlfNGet(1, mlfGca(NULL), _mxarray32_, NULL));
    /*
     * 
     * cyclefr=[];
     */
    mlfAssign(&cyclefr, _mxarray34_);
    /*
     * for i=1:length(lines)
     */
    {
        int v_ = mclForIntStart(1);
        int e_ = mclLengthInt(mclVv(lines, "lines"));
        if (v_ > e_) {
            mlfAssign(&i, _mxarray34_);
        } else {
            /*
             * if (get(lines(i),'ZData'))
             * xd=get(lines(i),'XData');
             * cyclefr=[cyclefr;xd(1)];
             * end
             * end
             */
            for (; ; ) {
                if (mlfTobool(
                      mlfNGet(
                        1,
                        mclIntArrayRef1(mclVv(lines, "lines"), v_),
                        _mxarray9_,
                        NULL))) {
                    mlfAssign(
                      &xd,
                      mlfNGet(
                        1,
                        mclIntArrayRef1(mclVv(lines, "lines"), v_),
                        _mxarray35_,
                        NULL));
                    mlfAssign(
                      &cyclefr,
                      mlfVertcat(
                        mclVv(cyclefr, "cyclefr"),
                        mclIntArrayRef1(mclVv(xd, "xd"), 1),
                        NULL));
                }
                if (v_ == e_) {
                    break;
                }
                ++v_;
            }
            mlfAssign(&i, mlfScalar(v_));
        }
    }
    /*
     * 
     * cyclefr=fix(cyclefr);
     */
    mlfAssign(&cyclefr, mlfFix(mclVv(cyclefr, "cyclefr")));
    /*
     * cyclefr(find(cyclefr<1))=[];
     */
    mlfIndexDelete(
      &cyclefr,
      "(?)",
      mlfFind(NULL, NULL, mclLt(mclVv(cyclefr, "cyclefr"), _mxarray2_)));
    /*
     * 
     * cyclefr=sort(cyclefr);
     */
    mlfAssign(&cyclefr, mlfSort(NULL, mclVv(cyclefr, "cyclefr"), NULL));
    mclValidateOutput(cyclefr, 1, nargout_, "cyclefr", "detectstepsplot");
    mclValidateOutput(*figh, 2, nargout_, "figh", "detectstepsplot");
    mxDestroyArray(fig);
    mxDestroyArray(ans);
    mxDestroyArray(cyclelines);
    mxDestroyArray(msg);
    mxDestroyArray(lines);
    mxDestroyArray(i);
    mxDestroyArray(xd);
    mxDestroyArray(mdf);
    mxDestroyArray(stpfr);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return cyclefr;
    /*
     * 
     * % uiwait(msgbox({'Beginning of cycles found at frame #',...
     * %                sprintf('%d  ',cyclefr)},'Cycles'))
     */
}

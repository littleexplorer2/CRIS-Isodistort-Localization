/* DERIVED GHIDRA PSEUDOCODE; NOT ORIGINAL SOURCE; NOT COMPILABLE CODE.
 * Evidence: ev_f52b50abf8a82a8fdff67c6627669397e75935686ac1fa0b0ba3e89d81c6b0e0
 * Raw response: evidence/native/rea-linux-iso-recovered/08-constant.json
 * Types, array bounds and control flow may be imperfectly recovered.
 */

undefined8 constant_(int *param_1)

{
  if (const_common_ == 0) {
    const_common_ = 1;
    init_constants_();
  }
  if ((*param_1 < 1) || (DAT_01bb1244 < *param_1)) {
    bomb2_("constants 52",0xc);
  }
  return *(undefined8 *)(&const_common_ + (long)*param_1 * 2);
}


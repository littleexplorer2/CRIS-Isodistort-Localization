/* DERIVED GHIDRA PSEUDOCODE; NOT ORIGINAL SOURCE; NOT COMPILABLE CODE.
 * Evidence: ev_532c54c72d30b144fda62477b3ffc34e3b611d5054f43e5e06a71ad360dd6eaa
 * Raw response: evidence/native/rea-linux-iso-recovered/05-orderparam-check.json
 * Types, array bounds and control flow may be imperfectly recovered.
 */

undefined8 orderparam_check_(int *param_1,int *param_2,long param_3,double *param_4)

{
  int iVar1;
  int iVar2;
  double *pdVar3;
  int iVar4;
  double *pdVar5;
  long lVar6;
  long lVar7;
  long lVar8;
  long local_60;
  int local_50;
  double local_40 [2];
  
  iVar1 = *param_2;
  if (0 < iVar1) {
    local_60 = 0x30;
    local_50 = 1;
    pdVar3 = param_4;
    do {
      iVar2 = *param_1;
      if (0 < iVar2) {
        lVar7 = 1;
        lVar8 = param_3;
        do {
          lVar8 = lVar8 + 8;
          if (*param_1 < 1) {
            local_40[0] = 0.0;
          }
          else {
            local_40[0] = 0.0;
            pdVar5 = pdVar3;
            lVar6 = lVar8;
            do {
              local_40[0] = local_40[0] + *(double *)(lVar6 + -8) * *pdVar5;
              lVar6 = lVar6 + 0x180;
              pdVar5 = pdVar5 + 1;
            } while (pdVar5 != param_4 + (ulong)(*param_1 - 1) + local_60 + -0x2f);
          }
          local_40[0] = local_40[0] - pdVar3[lVar7 + -1];
          iVar4 = ncmp_(local_40);
          if (iVar4 != 0) {
            return 0;
          }
          lVar7 = lVar7 + 1;
        } while ((int)lVar7 <= iVar2);
      }
      local_50 = local_50 + 1;
      local_60 = local_60 + 0x30;
      pdVar3 = pdVar3 + 0x30;
    } while (local_50 <= iVar1);
  }
  return 1;
}


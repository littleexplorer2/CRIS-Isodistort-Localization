/* DERIVED GHIDRA PSEUDOCODE; NOT ORIGINAL SOURCE; NOT COMPILABLE CODE.
 * Evidence: ev_995c4e8a0963d7a321775c56b99b779feec99547d7d3799555092da781337c30
 * Raw response: evidence/native/rea-linux-iso-recovered/13-project-vector.json
 * Types, array bounds and control flow may be imperfectly recovered.
 */

void project_vector_(int *param_1,int *param_2,int *param_3,int *param_4,int *param_5,int *param_6,
                    int *param_7,int *param_8,long param_9,int *param_10,long param_11)

{
  double dVar1;
  int iVar2;
  int iVar3;
  int iVar4;
  double *pdVar5;
  undefined4 *puVar6;
  long lVar7;
  long lVar8;
  long lVar9;
  long lVar10;
  double dVar11;
  double dVar12;
  double local_d0;
  long local_c0;
  int local_b0;
  int local_98;
  int local_94;
  long local_90;
  long local_88;
  int local_80;
  double local_60 [4];
  int local_3c [3];
  
  local_d0 = 1.0;
  if ((*param_6 == 2) && (*(int *)(&DAT_00e7121c + (long)*param_7 * 4) == 0)) {
    local_d0 = -1.0;
  }
  local_94 = 0;
  local_98 = 1;
  do {
    local_80 = local_94;
    if (*param_3 ==
        (int)*(char *)((long)local_98 + (long)*param_6 * 3 + (long)*param_1 * 6 + 0x1a04e1d)) {
      iVar4 = *param_8;
      if (iVar4 == 0) {
        local_80 = local_94 + 1;
      }
      else if (0 < iVar4) {
        local_90 = param_9 + 0x498;
        local_88 = ((long)local_94 * 9 + 9) * 0x10;
        local_60[0] = (double)(param_11 + 0x18);
        do {
          local_80 = local_80 + 1;
          iVar3 = *param_4;
          if (0 < iVar3) {
            local_c0 = (long)local_60[0] + local_88 * 8;
            local_b0 = 1;
            lVar10 = local_90;
            do {
              local_60[1] = 0.0;
              local_60[2] = 0.0;
              local_60[3] = 0.0;
              iVar2 = *param_5;
              if (0 < iVar2) {
                lVar8 = 1;
                lVar9 = (long)local_98 * 3;
                do {
                  lVar7 = 1;
                  do {
                    local_3c[0] = (int)*(char *)(lVar7 + (long)*param_6 * 9 + (long)*param_2 * 0x12
                                                         + 0x1a04ec8 + lVar9);
                    dVar12 = local_60[lVar7];
                    dVar1 = *(double *)(lVar10 + -0x4a0 + lVar8 * 8);
                    dVar11 = (double)constant_(local_3c);
                    local_60[lVar7] = dVar11 * dVar1 + dVar12;
                    lVar7 = lVar7 + 1;
                  } while (lVar7 != 4);
                  lVar8 = lVar8 + 1;
                  lVar9 = lVar9 + 3;
                } while ((int)lVar8 <= iVar2);
              }
              iVar2 = *param_7;
              lVar9 = 1;
              do {
                puVar6 = &common_space_ + (long)iVar2 * 9 + lVar9;
                pdVar5 = local_60 + 1;
                dVar12 = 0.0;
                do {
                  dVar12 = dVar12 + (double)(int)puVar6[-10] * local_d0 * *pdVar5;
                  puVar6 = puVar6 + 3;
                  pdVar5 = pdVar5 + 1;
                } while (pdVar5 != (double *)&stack0xffffffffffffffc0);
                *(double *)(local_c0 + -0x4a0 + lVar9 * 8) = dVar12;
                lVar9 = lVar9 + 1;
              } while (lVar9 != 4);
              local_b0 = local_b0 + 1;
              local_c0 = local_c0 + 0x18;
              lVar10 = lVar10 + 0x18;
            } while (local_b0 <= iVar3);
          }
          local_90 = local_90 + 0x480;
          local_88 = local_88 + 0x90;
        } while ((local_80 - local_94) + 1 <= iVar4);
      }
    }
    local_94 = local_80;
    local_98 = local_98 + *param_5;
  } while (local_98 < 4);
  *param_10 = local_94;
  return;
}


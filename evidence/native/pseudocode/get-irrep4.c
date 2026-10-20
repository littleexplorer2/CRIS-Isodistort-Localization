/* DERIVED GHIDRA PSEUDOCODE; NOT ORIGINAL SOURCE; NOT COMPILABLE CODE.
 * Evidence: ev_14bd60fa1fca5963680c7f1cfd61949e8da31466cb5b6fab89b10b6f87edab71
 * Raw response: evidence/native/rea-linux-iso-recovered/10-get-irrep4.json
 * Types, array bounds and control flow may be imperfectly recovered.
 */

/* WARNING: Globals starting with '_' overlap smaller symbols at the same address */

void get_irrep4_(int *param_1,int *param_2,undefined8 *param_3,undefined8 param_4,double *param_5,
                void *param_6,int *param_7)

{
  char cVar1;
  short sVar2;
  short sVar3;
  short sVar4;
  int iVar5;
  undefined4 uVar6;
  undefined4 *puVar7;
  long lVar8;
  double *pdVar9;
  undefined1 *puVar10;
  long lVar11;
  int *piVar12;
  undefined4 *puVar13;
  int iVar14;
  undefined *puVar15;
  long lVar16;
  long lVar17;
  undefined8 *puVar18;
  long lVar19;
  undefined *puVar20;
  byte bVar21;
  double dVar22;
  double local_49a8 [49];
  double adStack_4820 [2255];
  double local_1a8 [24];
  undefined1 local_e8 [12];
  int local_dc;
  undefined4 local_d8;
  undefined4 local_d4;
  undefined4 local_d0;
  int local_cc;
  undefined8 local_c8 [12];
  int local_68;
  undefined8 local_58;
  undefined8 uStack_50;
  int local_48;
  
  bVar21 = 0;
  local_58 = *param_3;
  uStack_50 = param_3[1];
  local_48 = *(int *)(param_3 + 2);
  iVar5 = *param_2;
  if (iVar5 != irrpointersave_9) {
    sVar2 = *(short *)(common_little_ + (long)iVar5 * 2 + 0x2d6f6e);
    ispace_8 = (int)sVar2;
    cVar1 = common_little_[(long)iVar5 + 0x36c333];
    kpnt_7 = (int)cVar1;
    latt_6 = (int)(char)(&DAT_00e74de3)[sVar2];
    lVar16 = (long)sVar2 * 0x1b + (long)cVar1;
    sVar3 = *(short *)(common_little_ + lVar16 * 2 + 0x2de502);
    nmod_5 = (int)sVar3;
    sVar4 = *(short *)(common_little_ + lVar16 * 2 + 0x2e1586);
    irrpointersave_9 = iVar5;
    if (*(char *)((long)cVar1 + (long)(char)(&DAT_00e74de3)[sVar2] * 0x1b + 0x119db56) == '\0') {
      generate_space_group_(&ispace_8,&nops_4,iops_3);
      puVar18 = &mapops_2;
      for (lVar16 = 0x24; lVar16 != 0; lVar16 = lVar16 + -1) {
        *puVar18 = 0;
        puVar18 = puVar18 + (ulong)bVar21 * -2 + 1;
      }
      if (0 < nops_4) {
        piVar12 = (int *)(iops_3 + 0x10);
        iVar14 = nops_4 + 1;
        iVar5 = 1;
        do {
          *(int *)(&DAT_00e0afdc + (long)*piVar12 * 4) = iVar5;
          iVar5 = iVar5 + 1;
          piVar12 = piVar12 + 5;
        } while (iVar5 != iVar14);
      }
      lVar16 = (long)sVar4;
      if (0 < nmod_5) {
        lVar19 = lVar16 * 4 + -5;
        puVar10 = common_little_ + lVar16 * 0x10;
        puVar13 = &DAT_00e0aa20;
        lVar17 = (long)nmod_5;
        do {
          lVar11 = 0;
          do {
            *(undefined4 *)((long)puVar13 + lVar11 + -0x40) =
                 *(undefined4 *)(puVar10 + lVar11 + 0x24c498);
            lVar11 = lVar11 + 4;
          } while (lVar11 != 0x10);
          lVar19 = lVar19 + 4;
          puVar10 = puVar10 + 0x10;
          puVar13 = puVar13 + 0x10;
        } while (lVar19 != (lVar16 + lVar17) * 4 + -5);
        puVar13 = &DAT_00e0aa20;
        lVar16 = -5;
        do {
          puVar7 = puVar13 + -0xc;
          do {
            *puVar7 = 0;
            puVar7[1] = 0;
            puVar7[2] = 0;
            puVar7[3] = 0;
            puVar7 = puVar7 + 4;
          } while (puVar7 != puVar13);
          lVar16 = lVar16 + 0x10;
          puVar13 = puVar13 + 0x10;
        } while (lVar16 != lVar17 * 0x10 + -5);
      }
    }
    else if (0 < sVar3) {
      lVar16 = (long)sVar4 * 0x40 + 0xea1a3c;
      puVar20 = &DAT_00e0a9dc;
      lVar19 = 0xf;
      do {
        lVar17 = lVar19 + -0x10;
        lVar11 = lVar16;
        puVar15 = puVar20;
        do {
          lVar8 = 0;
          do {
            *(undefined4 *)(puVar15 + lVar8 + 4) = *(undefined4 *)(lVar11 + 0x5223c + lVar8);
            lVar8 = lVar8 + 4;
          } while (lVar8 != 0x10);
          lVar17 = lVar17 + 4;
          lVar11 = lVar11 + 0x10;
          puVar15 = puVar15 + 0x10;
        } while (lVar17 != lVar19);
        lVar19 = lVar19 + 0x10;
        lVar16 = lVar16 + 0x40;
        puVar20 = puVar20 + 0x40;
      } while (lVar19 != (long)sVar3 * 0x10 + 0xf);
    }
  }
  if (*(char *)((long)kpnt_7 + (long)latt_6 * 0x1b + 0x119db56) == '\0') {
    iVar5 = *(int *)(&DAT_00e0afdc + (long)local_48 * 4);
    if (iVar5 == 0) {
      bomb2_("get_irrep4 63",0xd);
    }
    vsub_(&local_58,iops_3 + ((long)iVar5 * 5 + -5) * 4,&local_d8);
  }
  else {
    local_d8 = (undefined4)local_58;
    local_d4 = local_58._4_4_;
    local_d0 = (undefined4)uStack_50;
    local_cc = uStack_50._4_4_;
  }
  puVar18 = &iopn_0;
  for (lVar16 = 0xc; iVar5 = nmod_5, lVar16 != 0; lVar16 = lVar16 + -1) {
    *puVar18 = 0;
    puVar18 = puVar18 + (ulong)bVar21 * -2 + 1;
  }
  _DAT_00e0a9c0 = 1;
  if (0 < nmod_5) {
    puVar13 = &kvec_1;
    lVar16 = 1;
    do {
      get_kvec_(puVar13,param_4,local_e8);
      puVar18 = local_c8;
      for (lVar19 = 0xc; lVar19 != 0; lVar19 = lVar19 + -1) {
        *puVar18 = 0;
        puVar18 = puVar18 + (ulong)bVar21 * -2 + 1;
      }
      uVar6 = ivdot_(&DAT_005b08e0,local_e8,&local_d8);
      *(undefined4 *)((long)local_c8 + lVar16 * 4 + -4) = uVar6;
      local_68 = local_dc * local_cc;
      vaddn_(local_c8,&iopn_0,&iopn_0,&DAT_005b08e4);
      lVar16 = lVar16 + 1;
      puVar13 = puVar13 + 0x10;
    } while ((int)lVar16 <= iVar5);
  }
  memset(local_49a8,0,0x4800);
  lVar16 = 0;
  do {
    local_1a8[lVar16] = (double)*(int *)((long)&iopn_0 + lVar16 * 4) / (double)_DAT_00e0a9c0;
    lVar16 = lVar16 + 1;
  } while (lVar16 != 0x18);
  get_irrep4_ssgn_(param_2,&local_58,local_1a8,local_49a8,param_7);
  if (*param_7 != 0) {
    bomb2_("get_irrep4 93",0xd);
  }
  if (*param_1 == 2) {
    memcpy(param_6,local_49a8,0x4800);
  }
  else if (*param_1 == 1) {
    if ((char)common_little_[(long)*param_2 + 0x369af7] < '\x01') {
      dVar22 = 0.0;
    }
    else {
      pdVar9 = local_49a8;
      dVar22 = 0.0;
      do {
        dVar22 = dVar22 + *pdVar9;
        pdVar9 = pdVar9 + 0x31;
      } while (adStack_4820 +
               (ulong)((int)(char)common_little_[(long)*param_2 + 0x369af7] - 1) * 0x31 != pdVar9);
    }
    *param_5 = dVar22;
  }
  else {
    bomb2_("get_irrep4 106",0xe);
  }
  return;
}


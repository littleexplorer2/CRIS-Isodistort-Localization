/* DERIVED GHIDRA PSEUDOCODE; NOT ORIGINAL SOURCE; NOT COMPILABLE CODE.
 * Evidence: ev_876f3330452ee37940cf889e2523bf929f49ca0cfd8c7a84ee870ada181f9d86
 * Raw response: evidence/native/rea-linux-iso-recovered/11-get-irreps.json
 * Types, array bounds and control flow may be imperfectly recovered.
 */

void get_irreps_(int *param_1,int *param_2,int *param_3,undefined8 *param_4,undefined8 param_5,
                double *param_6,long param_7,int *param_8,int *param_9)

{
  char cVar1;
  int iVar2;
  int iVar3;
  long lVar4;
  long lVar5;
  undefined1 *puVar6;
  double unaff_RBX;
  long lVar7;
  long lVar8;
  bool bVar9;
  undefined1 local_4868 [384];
  undefined1 local_46e8 [384];
  undefined1 auStack_4568 [17672];
  int local_60;
  int local_5c;
  undefined8 local_58;
  undefined8 uStack_50;
  int local_48;
  double local_40 [2];
  
  iVar2 = *param_8;
  iVar3 = *param_1;
  if (iVar3 - 1U < 2) {
    local_58 = *param_4;
    uStack_50 = param_4[1];
    local_48 = *(int *)(param_4 + 2);
    unaff_RBX = 1.0;
    local_60 = iVar3;
  }
  else if (iVar3 - 3U < 2) {
    local_58 = *param_4;
    uStack_50 = param_4[1];
    local_48 = (int)(char)common_magnetic_[(long)*(int *)(param_4 + 2) + 0x102078];
    unaff_RBX = -1.0;
    local_60 = iVar3 + -2;
    if (*(int *)(common_magnetic_ + (long)*(int *)(param_4 + 2) * 4 + 0x24fb4) == 0) {
      unaff_RBX = 1.0;
    }
  }
  else {
    bomb2_("get_irreps 74",0xd);
  }
  local_5c = *param_2;
  if (local_5c == 0) {
    iVar3 = *param_3;
    if (iVar3 == 0) {
      *param_9 = 1;
      return;
    }
    if (0x12a8 < iVar3 - 1U) {
      *param_9 = 1;
      return;
    }
    local_5c = (int)*(short *)(common_little_ + (long)iVar3 * 2 + 0x2dbfe6);
  }
  else if ((local_5c < 1) || ((short)common_little_._3069902_2_ < local_5c)) {
    *param_9 = 1;
    return;
  }
  lVar8 = (long)iVar2;
  if (lVar8 < 0) {
    lVar8 = 0;
  }
  cVar1 = common_little_[(long)local_5c + 0x369af7];
  if (local_60 == 2) {
    iVar2 = *param_8;
    if (iVar2 < cVar1) {
      *param_9 = 1;
      return;
    }
    if (0 < iVar2) {
      lVar5 = param_7 + -8;
      lVar7 = 0;
      do {
        lVar4 = 1;
        do {
          *(undefined8 *)(lVar5 + lVar4 * 8) = 0;
          bVar9 = lVar4 != iVar2;
          lVar4 = lVar4 + 1;
        } while (bVar9);
        lVar7 = lVar7 + 1;
        lVar5 = lVar5 + lVar8 * 8;
      } while (lVar7 != iVar2);
    }
  }
  get_irrep4_(&local_60,&local_5c,&local_58,param_5,local_40,local_4868,param_9);
  if (*param_9 == 0) {
    if (local_60 == 1) {
      *param_6 = unaff_RBX * local_40[0];
    }
    else if ('\0' < cVar1) {
      puVar6 = local_46e8;
      param_7 = param_7 + -8;
      do {
        lVar5 = 1;
        do {
          *(double *)(param_7 + lVar5 * 8) = unaff_RBX * *(double *)(puVar6 + lVar5 * 8 + -0x188);
          lVar5 = lVar5 + 1;
        } while ((int)lVar5 <= (int)cVar1);
        puVar6 = puVar6 + 0x180;
        param_7 = param_7 + lVar8 * 8;
      } while (puVar6 != auStack_4568 + (ulong)((int)cVar1 - 1) * 0x180);
    }
  }
  return;
}


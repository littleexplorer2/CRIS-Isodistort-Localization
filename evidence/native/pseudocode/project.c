/* DERIVED GHIDRA PSEUDOCODE; NOT ORIGINAL SOURCE; NOT COMPILABLE CODE.
 * Evidence: ev_1131e1fb43cc1bbde64db04c22132ae5b6ee430b157e4ce0589a9f3a02ca8cc5
 * Raw response: evidence/native/rea-linux-iso-recovered/12-project.json
 * Types, array bounds and control flow may be imperfectly recovered.
 */

/* WARNING: Type propagation algorithm not settling */

void project_(int *param_1,undefined8 param_2,int *param_3,long param_4,int *param_5,int *param_6,
             int *param_7,long param_8)

{
  char cVar1;
  char cVar2;
  char cVar3;
  char cVar4;
  short sVar5;
  short sVar6;
  bool bVar7;
  int iVar8;
  int iVar9;
  int iVar10;
  int iVar11;
  int iVar12;
  int *piVar13;
  ulong uVar14;
  long lVar15;
  double *pdVar16;
  double *pdVar17;
  long lVar18;
  double *pdVar19;
  bool bVar20;
  int iVar21;
  long lVar22;
  long lVar23;
  long lVar24;
  long lVar25;
  int iVar26;
  int iVar27;
  long lVar28;
  double *pdVar29;
  double *pdVar30;
  int iVar31;
  double *pdVar32;
  undefined1 *puVar33;
  long lVar34;
  int iVar35;
  double dVar36;
  double dVar37;
  int iStack_17ebc;
  long lStack_17eb8;
  double dStack_17e88;
  undefined1 auStack_17e38 [36096];
  double adStack_f138 [96];
  double local_ee38 [4609];
  double local_5e30;
  undefined1 local_5e28 [18424];
  double adStack_1630 [578];
  double local_420;
  int local_418;
  int local_414;
  undefined1 local_410 [4];
  int local_40c;
  int local_408 [5];
  int aiStack_3f4 [238];
  int local_3c [3];
  
  lVar22 = (long)*param_1;
  sVar5 = *(short *)(common_little_ + lVar22 * 2 + 0x2d1ef6);
  sVar6 = *(short *)(common_little_ + lVar22 * 2 + 0x2d6f6e);
  cVar1 = common_little_[lVar22 + 0x369af7];
  cVar2 = common_little_[lVar22 + 0x35f9e7];
  if (cVar2 != '\x01') {
    transform_irrep2complex_(param_1,auStack_17e38,local_ee38,&DAT_005afb80);
  }
  iVar8 = *param_3;
  lVar22 = (long)(char)common_wyckoff_
                       [(long)(iVar8 + -1 +
                              (int)*(short *)(common_wyckoff_ + (long)sVar6 * 2 + 0x2fa22)) +
                        0x559fb];
  iVar35 = *param_5;
  local_40c = (int)*(short *)(common_wyckoff_ + ((long)iVar35 + 0x274c6 + lVar22 * 0xc) * 2);
  cVar3 = common_image_
          [(long)*(short *)(common_irrep_ +
                           (long)*(short *)(common_wyckoff_ +
                                           ((long)iVar35 + 0x274c6 + lVar22 * 0xc) * 2) * 2 +
                           0x551f4) + 0x20e71];
  if (sVar5 == 0) {
    lVar23 = (long)(iVar8 + *(int *)(common_little_ + (long)*param_1 * 4 + 0x6104) + -1);
    iVar8 = *(int *)(common_little_ + lVar23 * 4 + 0x101f4);
    if (iVar8 < (char)common_little_[lVar23 + 0x3713ab] + iVar8) {
      lVar24 = (long)iVar8 + -1;
      do {
        if (iVar35 == (char)common_little_[lVar24 + 0x39656c]) {
          cVar4 = common_little_[lVar24 + 0x37cb44];
          goto LAB_0049f6a3;
        }
        lVar24 = lVar24 + 1;
      } while (lVar24 != (long)iVar8 + (ulong)((int)(char)common_little_[lVar23 + 0x3713ab] - 1));
    }
    *param_7 = 0;
  }
  else {
    lVar23 = (long)(iVar8 + *(int *)(&DAT_0199367c + (long)sVar5 * 4) + -1);
    iVar8 = *(int *)(common_wyckoff_ + lVar23 * 4 + 0x4aa0);
    if (iVar8 < (char)common_wyckoff_[lVar23 + 0x57d5e] + iVar8) {
      lVar24 = (long)iVar8 + -1;
      do {
        if (iVar35 == *(short *)(common_wyckoff_ + lVar24 * 2 + 0x30d16)) {
          cVar4 = common_wyckoff_[lVar24 + 0x6293f];
LAB_0049f6a3:
          iVar8 = (int)cVar4;
          if (common_image_
              [(long)*(short *)(common_irrep_ +
                               (long)*(short *)(common_wyckoff_ +
                                               ((long)iVar35 + 0x274c6 + lVar22 * 0xc) * 2) * 2 +
                               0x551f4) + 0x21587] != '\x01') {
            iVar8 = cVar4 * 2;
          }
          iVar35 = *(int *)(&DAT_005afb94 + (long)cVar2 * 4);
          if (iVar8 % iVar35 != 0) {
            bomb2_("project 149",0xb);
          }
          iVar8 = iVar8 / iVar35;
          *param_7 = iVar8;
          cVar4 = (&DAT_00e748e9)[lVar22];
          if (cVar4 < '\x01') {
            iVar31 = (int)*(short *)(common_wyckoff_ + (lVar22 + -1) * 2 + 0x4eca6);
            iVar35 = iVar31;
            if ((char)common_wyckoff_[lVar22 + 0x71786] < '\x01') goto LAB_0049f879;
          }
          else {
            piVar13 = local_408;
            uVar14 = (ulong)((int)cVar4 - 1);
            do {
              *piVar13 = 0;
              piVar13[1] = 0;
              piVar13[2] = 0;
              piVar13[3] = 1;
              piVar13 = piVar13 + 5;
            } while (piVar13 != aiStack_3f4 + uVar14 * 5);
            iVar31 = (int)*(short *)(common_wyckoff_ + (lVar22 + -1) * 2 + 0x4eca6);
            iVar26 = (char)common_wyckoff_[lVar22 + 0x71786] + iVar31;
            iVar35 = iVar31;
            if (iVar31 < iVar26) {
              lVar23 = (long)*(short *)(common_wyckoff_ + (lVar22 + -1) * 2 + 0x4eca6) * 0x30;
              lVar22 = lVar23;
LAB_0049f84e:
              puVar33 = common_wyckoff_ + lVar22;
              do {
                lVar24 = 4;
                while (*(int *)(param_4 + lVar24 * 4) != (int)(char)puVar33[0x5608f]) {
                  lVar24 = lVar24 + 5;
                  if (lVar24 == uVar14 * 5 + 9) {
                    iVar31 = iVar31 + 1;
                    lVar22 = lVar22 + 0x30;
                    iVar35 = iVar26;
                    if (iVar26 != iVar31) goto LAB_0049f84e;
                    goto LAB_0049f879;
                  }
                }
                local_408[lVar24] = (int)(char)puVar33[(lVar23 + 0x5608f) - lVar22];
                puVar33 = puVar33 + 1;
              } while (puVar33 != (undefined1 *)(uVar14 + 0x1993681 + lVar22));
            }
            else {
LAB_0049f879:
              iVar31 = iVar35;
              bomb2_("project 190",0xb);
            }
          }
          iVar26 = (int)cVar1;
          iVar35 = (int)cVar3;
          *param_6 = iVar31;
          local_418 = 0;
          if (cVar2 == '\x01') {
            local_414 = iVar26 * iVar35;
          }
          else {
            local_414 = (iVar26 * iVar35) / 2;
          }
          if (cVar1 < '\x01') goto LAB_004a0155;
          iVar9 = iVar26 / 2;
          iVar31 = iVar9 + 1;
          uVar14 = (ulong)((int)cVar1 - 1);
          lStack_17eb8 = 0x30;
          iStack_17ebc = 0;
          iVar10 = (int)cVar3;
          do {
            lVar23 = (long)(iStack_17ebc + 1);
            lVar24 = lVar23 * 0x90;
            lVar22 = lVar24 + 3;
            pdVar19 = (double *)(param_8 + -0x498 + lVar22 * 8);
            pdVar32 = (double *)(param_8 + -0x490 + ((ulong)(iVar10 - 1) + 3 + lVar24) * 8);
            pdVar29 = (double *)(param_8 + -0x468 + (lVar24 + uVar14 * 3) * 8);
            pdVar30 = pdVar19;
            pdVar17 = pdVar32;
            do {
              pdVar16 = pdVar30;
              if ('\0' < cVar3) {
                do {
                  *pdVar16 = 0.0;
                  pdVar16 = pdVar16 + 1;
                } while (pdVar16 != pdVar17);
              }
              pdVar30 = pdVar30 + 3;
              pdVar17 = pdVar17 + 3;
            } while (pdVar30 != pdVar29);
            if ('\0' < cVar4) {
              piVar13 = local_408;
              lVar28 = lVar22 * 8 + param_8;
              lVar34 = param_4;
              do {
                get_irreps_(&DAT_005afb88,param_1,&DAT_005afb84,lVar34,param_2,&local_420,local_5e28
                            ,&DAT_005afb80,local_3c);
                if (local_3c[0] != 0) {
                  bomb2_("project 223",0xb);
                }
                get_irreps_(&DAT_005afb88,&DAT_005afb84,&local_40c,piVar13,local_410,&local_420,
                            adStack_1630 + 1,&DAT_005afb8c,local_3c);
                lVar18 = lVar22;
                lVar25 = lVar28;
                puVar33 = local_5e28 + lStack_17eb8 * 8;
                if (local_3c[0] != 0) {
                  bomb2_("project 227",0xb);
                }
                do {
                  if ('\0' < cVar3) {
                    dVar37 = *(double *)(puVar33 + -0x180);
                    lVar15 = 1;
                    do {
                      *(double *)(lVar25 + -0x4a0 + lVar15 * 8) =
                           dVar37 * adStack_1630[lVar15] + *(double *)(lVar25 + -0x4a0 + lVar15 * 8)
                      ;
                      lVar15 = lVar15 + 1;
                    } while ((int)lVar15 <= iVar35);
                  }
                  lVar18 = lVar18 + 3;
                  lVar25 = lVar25 + 0x18;
                  puVar33 = puVar33 + 8;
                } while (lVar18 != uVar14 * 3 + 6 + lVar24);
                lVar34 = lVar34 + 0x14;
                piVar13 = piVar13 + 5;
              } while (lVar34 != param_4 + 0x14 + (ulong)((int)cVar4 - 1) * 0x14);
            }
            local_5e30 = 0.0;
            bVar7 = false;
            bVar20 = false;
            pdVar30 = pdVar32;
            pdVar17 = pdVar19;
            dVar37 = 0.0;
            do {
              pdVar16 = pdVar17;
              dVar36 = dVar37;
              if ('\0' < cVar3) {
                do {
                  dStack_17e88 = ABS(*pdVar16);
                  dVar36 = dStack_17e88;
                  if (dStack_17e88 <= dVar37) {
                    dVar36 = dVar37;
                  }
                  bVar7 = true;
                  if (dVar37 < dStack_17e88) {
                    bVar20 = true;
                  }
                  pdVar16 = pdVar16 + 1;
                  dVar37 = dVar36;
                } while (pdVar16 != pdVar30);
              }
              pdVar30 = pdVar30 + 3;
              pdVar17 = pdVar17 + 3;
              dVar37 = dVar36;
            } while (pdVar17 != pdVar29);
            if (bVar7) {
              local_420 = dStack_17e88;
            }
            if (bVar20) {
              local_5e30 = dVar36;
            }
            iVar11 = ncmp_(&local_5e30);
            pdVar30 = pdVar19;
            pdVar17 = pdVar32;
            if (iVar11 != 0) {
              do {
                pdVar16 = pdVar30;
                if ('\0' < cVar3) {
                  do {
                    *pdVar16 = *pdVar16 / local_5e30;
                    pdVar16 = pdVar16 + 1;
                  } while (pdVar16 != pdVar17);
                }
                pdVar30 = pdVar30 + 3;
                pdVar17 = pdVar17 + 3;
              } while (pdVar30 != pdVar29);
              if (cVar2 == '\x01') {
                iVar11 = 0;
                do {
                  if ('\0' < cVar3) {
                    pdVar17 = (double *)(zf_0 + ((long)local_418 + (long)iVar11 * 0x32) * 0x10);
                    pdVar30 = pdVar19;
                    do {
                      *pdVar17 = *pdVar30;
                      pdVar17[1] = 0.0;
                      pdVar30 = pdVar30 + 1;
                      pdVar17 = pdVar17 + 100;
                    } while (pdVar30 != pdVar32);
                    iVar11 = iVar11 + iVar10;
                  }
                  pdVar19 = pdVar19 + 3;
                  pdVar32 = pdVar32 + 3;
                  iVar27 = local_418 + 1;
                } while (pdVar19 != pdVar29);
LAB_0049fd15:
                local_418 = iVar27;
                zrowop2_(zf_0,&local_418,&local_414,&DAT_005afb94,&DAT_005afb90);
                if (0 < local_414) {
LAB_0049fd4c:
                  iVar11 = local_414;
                  lVar22 = 0x32;
                  iVar27 = 1;
                  do {
                    iVar12 = ncmpz_(zf_0 + (lVar22 + -0x33 + (long)local_418) * 0x10);
                    if (iVar12 != 0) {
                      iStack_17ebc = iStack_17ebc + 1;
                      if (iVar8 == iStack_17ebc) {
                        return;
                      }
                      goto LAB_0049fc20;
                    }
                    iVar27 = iVar27 + 1;
                    lVar22 = lVar22 + 0x32;
                  } while (iVar27 <= iVar11);
                  if (cVar2 == '\x02') goto LAB_004a01b3;
                }
                local_418 = local_418 + -1;
              }
              else {
                if (cVar2 != '\x02') {
                  iVar27 = local_418;
                  if ((cVar2 == '\x03') && (iVar27 = local_418 + 1, 1 < iVar26)) {
                    iVar12 = 0;
                    iVar11 = 1;
                    pdVar19 = adStack_f138 + (uVar14 + 0x30) * 2 + 2;
                    pdVar32 = local_ee38;
                    do {
                      if ('\0' < cVar3) {
                        pdVar30 = (double *)(zf_0 + ((long)local_418 + (long)iVar12 * 0x32) * 0x10);
                        iVar21 = 1;
                        pdVar29 = (double *)(param_8 + -0x480 + lVar23 * 0x480);
                        do {
                          dVar36 = 0.0;
                          dVar37 = 0.0;
                          pdVar17 = pdVar32;
                          pdVar16 = pdVar29;
                          do {
                            dVar37 = dVar37 + (*pdVar16 * *pdVar17 - pdVar17[1] * 0.0);
                            dVar36 = dVar36 + *pdVar17 * 0.0 + *pdVar16 * pdVar17[1];
                            pdVar16 = pdVar16 + 3;
                            pdVar17 = pdVar17 + 2;
                          } while (pdVar17 != pdVar19);
                          *pdVar30 = dVar37;
                          pdVar30[1] = dVar36;
                          pdVar30 = pdVar30 + 100;
                          pdVar29 = pdVar29 + 1;
                          iVar21 = iVar21 + 1;
                        } while (iVar21 <= iVar35);
                        iVar12 = iVar12 + iVar10;
                      }
                      iVar11 = iVar11 + 1;
                      pdVar32 = pdVar32 + 0x60;
                      pdVar19 = pdVar19 + 0x60;
                    } while (iVar11 <= iVar9);
                  }
                  goto LAB_0049fd15;
                }
                if (1 < iVar26) {
                  iVar11 = 0;
                  iVar27 = 1;
                  pdVar19 = adStack_f138 + (uVar14 + 0x30) * 2 + 2;
                  pdVar32 = local_ee38;
                  do {
                    if ('\0' < cVar3) {
                      pdVar30 = (double *)(zf_0 + ((long)local_418 + (long)iVar11 * 0x32) * 0x10);
                      iVar12 = 1;
                      pdVar29 = (double *)(param_8 + -0x480 + lVar23 * 0x480);
                      do {
                        dVar36 = 0.0;
                        dVar37 = 0.0;
                        pdVar17 = pdVar32;
                        pdVar16 = pdVar29;
                        do {
                          dVar37 = dVar37 + (*pdVar16 * *pdVar17 - pdVar17[1] * 0.0);
                          dVar36 = dVar36 + *pdVar17 * 0.0 + *pdVar16 * pdVar17[1];
                          pdVar16 = pdVar16 + 3;
                          pdVar17 = pdVar17 + 2;
                        } while (pdVar17 != pdVar19);
                        *pdVar30 = dVar37;
                        pdVar30[1] = dVar36;
                        pdVar30 = pdVar30 + 100;
                        pdVar29 = pdVar29 + 1;
                        iVar12 = iVar12 + 1;
                      } while (iVar12 <= iVar35);
                      iVar11 = iVar11 + iVar10;
                    }
                    iVar27 = iVar27 + 1;
                    pdVar32 = pdVar32 + 0x60;
                    pdVar19 = pdVar19 + 0x60;
                  } while (iVar27 <= iVar9);
                }
                if (iVar31 <= iVar26) {
                  iVar27 = 0;
                  pdVar19 = adStack_f138 + (long)iVar31 * 0x60;
                  iVar11 = iVar31;
                  do {
                    if ('\0' < cVar3) {
                      pdVar29 = (double *)
                                (zf_0 + ((long)local_418 + 1 + (long)iVar27 * 0x32) * 0x10);
                      iVar12 = 1;
                      pdVar32 = (double *)(param_8 + -0x480 + lVar23 * 0x480);
                      do {
                        dVar36 = 0.0;
                        dVar37 = 0.0;
                        iVar21 = 1;
                        pdVar30 = pdVar19;
                        pdVar17 = pdVar32;
                        do {
                          dVar37 = dVar37 + (*pdVar17 * *pdVar30 - pdVar30[1] * 0.0);
                          dVar36 = dVar36 + *pdVar30 * 0.0 + *pdVar17 * pdVar30[1];
                          iVar21 = iVar21 + 1;
                          pdVar17 = pdVar17 + 3;
                          pdVar30 = pdVar30 + 2;
                        } while (iVar21 <= iVar26);
                        *pdVar29 = dVar37;
                        pdVar29[1] = dVar36;
                        pdVar29 = pdVar29 + 100;
                        pdVar32 = pdVar32 + 1;
                        iVar12 = iVar12 + 1;
                      } while (iVar12 <= iVar35);
                      iVar27 = iVar27 + iVar10;
                    }
                    iVar11 = iVar11 + 1;
                    pdVar19 = pdVar19 + 0x60;
                  } while (iVar11 <= iVar26);
                }
                local_418 = local_418 + 2;
                zrowop2_(zf_0,&local_418,&local_414,&DAT_005afb94,&DAT_005afb90);
                if (0 < local_414) goto LAB_0049fd4c;
LAB_004a01b3:
                local_418 = local_418 + -2;
              }
            }
LAB_0049fc20:
            lStack_17eb8 = lStack_17eb8 + 0x30;
            if (lStack_17eb8 == uVar14 * 0x30 + 0x60) {
LAB_004a0155:
              bomb2_("project 353",0xb);
              return;
            }
          } while( true );
        }
        lVar24 = lVar24 + 1;
      } while (lVar24 != (long)iVar8 + (ulong)((int)(char)common_wyckoff_[lVar23 + 0x57d5e] - 1));
    }
    *param_7 = 0;
  }
  return;
}


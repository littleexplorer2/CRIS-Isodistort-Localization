/* DERIVED GHIDRA PSEUDOCODE; NOT ORIGINAL SOURCE; NOT COMPILABLE CODE.
 * Evidence: ev_2310ee83c419aac9e2bae0bac220527bd066f4f3a5ac69c79e2eca5861c873b7
 * Raw response: evidence/native/rea-linux-iso-recovered/14-get-isotropy-file.json
 * Types, array bounds and control flow may be imperfectly recovered.
 */

void get_isotropy_file_(int *param_1,int *param_2,undefined8 param_3,undefined8 *param_4,
                       void *param_5,undefined4 *param_6,int *param_7,undefined8 param_8,
                       undefined8 param_9,size_t param_10,size_t param_11)

{
  char cVar1;
  int iVar2;
  long lVar3;
  long lVar4;
  undefined4 local_268;
  undefined4 local_264;
  char *local_260;
  undefined4 local_258;
  undefined8 local_220;
  char *local_218;
  undefined8 local_210;
  undefined1 *local_1f8;
  undefined8 local_1f0;
  int local_58;
  int local_54;
  undefined1 local_50 [4];
  undefined1 local_4c [4];
  int local_48;
  int local_44;
  int local_40;
  undefined1 local_39 [9];
  
  iVar2 = *param_1;
  if (20000 < iVar2) {
    iVar2 = iVar2 + -20000;
  }
  *param_6 = 0;
  if (__isotropy_module_MOD_superspace == 0) {
    lVar3 = (long)iVar2;
    lVar4 = lVar3 + -1;
    local_44 = (int)*(short *)(common_little_ + lVar3 * 2 + 0x2d6f6e);
    cVar1 = *(char *)((long)(char)common_little_[lVar3 + 0x36c333] +
                     (long)(char)(&DAT_00e74de3)[*(short *)(common_little_ + lVar3 * 2 + 0x2d6f6e)]
                     * 0x1b + 0x119db56);
    local_48 = (int)cVar1;
    local_54 = (int)(char)common_little_[lVar3 + 0x369af7];
    iVar2 = *param_2;
    if (iVar2 == 1) {
      if (cVar1 == '\0') {
        if (0 < (long)param_10) {
          memset(param_4,0x20,param_10);
        }
      }
      else {
        if (0 < (long)param_10) {
          if ((long)param_10 < 0xd) {
            memcpy(param_4,"i     xx.iso",param_10);
          }
          else {
            *param_4 = 0x7878202020202069;
            *(undefined4 *)(param_4 + 1) = 0x6f73692e;
            memset((undefined1 *)((long)param_4 + 0xc),0x20,param_10 - 0xc);
          }
        }
        local_260 = "../source/get_isotropy_file.f";
        local_258 = 0x3e;
        local_1f8 = (undefined1 *)((long)param_4 + 1);
        local_1f0 = 5;
        local_220 = 0;
        local_264 = 0xffffffff;
        local_218 = "(i5.5)";
        local_210 = 6;
        local_268 = 0x5000;
        _gfortran_st_write(&local_268);
        _gfortran_transfer_integer_write(&local_268,param_1,4);
        _gfortran_st_write_done(&local_268);
        get_isotropy_filename_
                  (&local_48,param_3,param_4,&local_40,&local_58,local_4c,local_4c,local_50,param_10
                  );
        if (local_58 == 0) {
          if (local_40 == 0) {
            if (__isotropy_module_MOD_www == 0) {
              local_260 = "../source/get_isotropy_file.f";
              local_258 = 0x4e;
              local_218 = "(1x,a)";
              local_210 = 6;
              local_268 = 0x1000;
              local_264 = 6;
              _gfortran_st_write(&local_268);
              _gfortran_transfer_character_write
                        (&local_268,"You have requested information about isotropy subgroups for: ",
                         0x3d);
              _gfortran_st_write_done(&local_268);
              local_260 = "../source/get_isotropy_file.f";
              local_258 = 0x51;
              local_218 = "(6x,a,a,a,a)";
              local_210 = 0xc;
              local_268 = 0x1000;
              local_264 = 6;
              _gfortran_st_write(&local_268);
              _gfortran_transfer_character_write(&local_268,"irrep: ",7);
              _gfortran_transfer_character_write(&local_268,lVar4 * 8 + 0x1509ae6,8);
              _gfortran_transfer_character_write(&local_268," space group: ",0xe);
              _gfortran_transfer_character_write(&local_268,(long)local_44 * 8 + 0xe976fe,8);
              _gfortran_st_write_done(&local_268);
              local_260 = "../source/get_isotropy_file.f";
              local_258 = 0x53;
              local_218 = "(1x,a)";
              local_210 = 6;
              local_268 = 0x1000;
              local_264 = 6;
              _gfortran_st_write(&local_268);
              _gfortran_transfer_character_write
                        (&local_268,"The data base for these isotropy subgroups cannot be found.",
                         0x3b);
              _gfortran_st_write_done(&local_268);
              local_260 = "../source/get_isotropy_file.f";
              local_258 = 0x54;
              local_218 = "(1x,a)";
              local_210 = 6;
              local_268 = 0x1000;
              local_264 = 6;
              _gfortran_st_write(&local_268);
              _gfortran_transfer_character_write(&local_268,"Should the data base be added?",0x1e);
              _gfortran_st_write_done(&local_268);
              local_260 = "../source/get_isotropy_file.f";
              local_258 = 0x56;
              local_218 = "(1x,a)";
              local_210 = 6;
              local_268 = 0x1000;
              local_264 = 6;
              _gfortran_st_write(&local_268);
              _gfortran_transfer_character_write
                        (&local_268,"Enter RETURN to continue.  Enter any character to stop",0x36);
              _gfortran_st_write_done(&local_268);
              readinput_(local_39,1);
              lVar4 = _gfortran_string_len_trim(1,local_39);
              if (lVar4 != 0) {
                local_260 = "../source/get_isotropy_file.f";
                local_258 = 0x59;
                local_218 = "(a)";
                local_210 = 3;
                local_268 = 0x1000;
                local_264 = 6;
                _gfortran_st_write(&local_268);
                _gfortran_transfer_character_write(&local_268,"Quit display",0xc);
                _gfortran_st_write_done(&local_268);
                *param_6 = 1;
                return;
              }
            }
            local_260 = "../source/get_isotropy_file.f";
            local_258 = 0x5f;
            local_218 = "(1x,a19)";
            local_210 = 8;
            local_268 = 0x1000;
            local_264 = 6;
            _gfortran_st_write(&local_268);
            _gfortran_transfer_character_write(&local_268,"Adding data base...",0x13);
            _gfortran_st_write_done(&local_268);
            get_isotropy_(&local_44,param_1,param_3,&DAT_005a6020,local_4c,local_50);
          }
        }
        else {
          local_260 = "../source/get_isotropy_file.f";
          local_258 = 0x43;
          local_218 = "(1x,a)";
          local_210 = 6;
          local_268 = 0x1000;
          local_264 = 6;
          _gfortran_st_write(&local_268);
          _gfortran_transfer_character_write
                    (&local_268,"You have requested information about isotropy subgroups for: ",0x3d
                    );
          _gfortran_st_write_done(&local_268);
          local_260 = "../source/get_isotropy_file.f";
          local_258 = 0x46;
          local_218 = "(6x,a7,a8,a14,a8)";
          local_210 = 0x11;
          local_268 = 0x1000;
          local_264 = 6;
          _gfortran_st_write(&local_268);
          _gfortran_transfer_character_write(&local_268,"irrep: ",7);
          _gfortran_transfer_character_write(&local_268,lVar4 * 8 + 0x1509ae6,8);
          _gfortran_transfer_character_write(&local_268," space group: ",0xe);
          _gfortran_transfer_character_write(&local_268,(long)local_44 * 8 + 0xe976fe,8);
          _gfortran_st_write_done(&local_268);
          local_260 = "../source/get_isotropy_file.f";
          local_258 = 0x48;
          local_218 = "(1x,a)";
          local_210 = 6;
          local_268 = 0x1000;
          local_264 = 6;
          _gfortran_st_write(&local_268);
          _gfortran_transfer_character_write
                    (&local_268,"This data base is already being generated.  Try again later.",0x3c)
          ;
          _gfortran_st_write_done(&local_268);
          *param_6 = 1;
        }
      }
    }
    else if (iVar2 == 2) {
      if (0 < (long)param_11) {
        memset(param_5,0x20,param_11);
      }
      orderparam_label_(&local_54,&DAT_005a6024,param_5,param_11);
      if (local_48 == 0) {
        if (0 < (long)param_10) {
          memset(param_4,0x20,param_10);
        }
      }
      else {
        if (0 < (long)param_10) {
          if ((long)param_10 < 0xd) {
            memcpy(param_4,"i     xx.iso",param_10);
          }
          else {
            *param_4 = 0x7878202020202069;
            *(undefined4 *)(param_4 + 1) = 0x6f73692e;
            memset((undefined1 *)((long)param_4 + 0xc),0x20,param_10 - 0xc);
          }
        }
        local_260 = "../source/get_isotropy_file.f";
        local_258 = 0x70;
        local_1f8 = (undefined1 *)((long)param_4 + 1);
        local_1f0 = 5;
        local_220 = 0;
        local_264 = 0xffffffff;
        local_218 = "(i5.5)";
        local_210 = 6;
        local_268 = 0x5000;
        _gfortran_st_write(&local_268);
        _gfortran_transfer_integer_write(&local_268,param_1,4);
        _gfortran_st_write_done(&local_268);
        get_isotropy_filename_
                  (&local_48,param_3,param_4,&local_40,&local_58,local_4c,local_4c,local_50,param_10
                  );
        if ((local_40 == 0) || (local_58 != 0)) {
          *(undefined1 *)param_4 = 0x6b;
          get_isotropy_filename_
                    (&local_48,param_3,param_4,&local_40,&local_58,local_4c,local_4c,local_50,
                     param_10);
          if (local_58 == 0) {
            if (local_40 == 0) {
              get_isotropy_(&local_44,param_1,param_3,&DAT_005a6028,local_4c,local_50);
            }
          }
          else {
            local_260 = "../source/get_isotropy_file.f";
            local_258 = 0x79;
            local_218 = "(1x,a)";
            local_210 = 6;
            local_268 = 0x1000;
            local_264 = 6;
            _gfortran_st_write(&local_268);
            _gfortran_transfer_character_write
                      (&local_268,"You have requested information about the kernel for: ",0x35);
            _gfortran_st_write_done(&local_268);
            local_260 = "../source/get_isotropy_file.f";
            local_258 = 0x7c;
            local_218 = "(6x,a7,a8,a14,a8)";
            local_210 = 0x11;
            local_268 = 0x1000;
            local_264 = 6;
            _gfortran_st_write(&local_268);
            _gfortran_transfer_character_write(&local_268,"irrep: ",7);
            _gfortran_transfer_character_write(&local_268,lVar4 * 8 + 0x1509ae6,8);
            _gfortran_transfer_character_write(&local_268," space group: ",0xe);
            _gfortran_transfer_character_write(&local_268,(long)local_44 * 8 + 0xe976fe,8);
            _gfortran_st_write_done(&local_268);
            local_260 = "../source/get_isotropy_file.f";
            local_258 = 0x7e;
            local_218 = "(1x,a)";
            local_210 = 6;
            local_268 = 0x1000;
            local_264 = 6;
            _gfortran_st_write(&local_268);
            _gfortran_transfer_character_write
                      (&local_268,"This data base is already being generated.  Try again later.",
                       0x3c);
            _gfortran_st_write_done(&local_268);
            *param_6 = 1;
          }
        }
      }
    }
    else if (iVar2 - 3U < 2) {
      if (0 < (long)param_10) {
        if ((long)param_10 < 0xd) {
          memcpy(param_4,&DAT_005a5e24,param_10);
        }
        else {
          *param_4 = 0x7878202020202061;
          *(undefined4 *)(param_4 + 1) = 0x6f73692e;
          memset((undefined1 *)((long)param_4 + 0xc),0x20,param_10 - 0xc);
        }
      }
      local_260 = "../source/get_isotropy_file.f";
      local_258 = 0x89;
      local_1f8 = (undefined1 *)((long)param_4 + 1);
      local_1f0 = 5;
      local_220 = 0;
      local_264 = 0xffffffff;
      local_218 = "(i5.5)";
      local_210 = 6;
      local_268 = 0x5000;
      _gfortran_st_write(&local_268);
      _gfortran_transfer_integer_write(&local_268,param_1,4);
      _gfortran_st_write_done(&local_268);
      get_isotropy_filename_
                (&local_48,param_3,param_4,&local_40,&local_58,local_4c,local_4c,local_50,param_10);
      if (local_58 == 0) {
        if (local_40 == 0) {
          get_isotropy_(&local_44,param_1,param_3,&DAT_005a6024,local_4c,local_50);
        }
      }
      else {
        local_260 = "../source/get_isotropy_file.f";
        local_258 = 0x8d;
        local_218 = "(1x,a)";
        local_210 = 6;
        local_268 = 0x1000;
        local_264 = 6;
        _gfortran_st_write(&local_268);
        _gfortran_transfer_character_write
                  (&local_268,"You have requested information about isotropy subgroups for: ",0x3d);
        _gfortran_st_write_done(&local_268);
        local_260 = "../source/get_isotropy_file.f";
        local_258 = 0x90;
        local_218 = "(6x,a7,a8,a14,a8)";
        local_210 = 0x11;
        local_268 = 0x1000;
        local_264 = 6;
        _gfortran_st_write(&local_268);
        _gfortran_transfer_character_write(&local_268,"irrep: ",7);
        _gfortran_transfer_character_write(&local_268,lVar4 * 8 + 0x1509ae6,8);
        _gfortran_transfer_character_write(&local_268," space group: ",0xe);
        _gfortran_transfer_character_write(&local_268,(long)local_44 * 8 + 0xe976fe,8);
        _gfortran_st_write_done(&local_268);
        local_260 = "../source/get_isotropy_file.f";
        local_258 = 0x92;
        local_218 = "(1x,a)";
        local_210 = 6;
        local_268 = 0x1000;
        local_264 = 6;
        _gfortran_st_write(&local_268);
        _gfortran_transfer_character_write
                  (&local_268,"This data base is already being generated.  Try again later.",0x3c);
        _gfortran_st_write_done(&local_268);
        *param_6 = 1;
      }
    }
    else if (iVar2 == 5) {
      if (local_54 == *param_7) {
        if (0 < (long)param_11) {
          memset(param_5,0x20,param_11);
        }
        orderparam_label_(param_8,&DAT_005a6024,param_5,param_11);
        if (0 < (long)param_10) {
          if ((long)param_10 < 0xd) {
            memcpy(param_4,"v     xx.iso",param_10);
          }
          else {
            *param_4 = 0x7878202020202076;
            *(undefined4 *)(param_4 + 1) = 0x6f73692e;
            memset((undefined1 *)((long)param_4 + 0xc),0x20,param_10 - 0xc);
          }
        }
        local_260 = "../source/get_isotropy_file.f";
        local_258 = 0xaa;
        local_1f8 = (undefined1 *)((long)param_4 + 1);
        local_1f0 = 5;
        local_220 = 0;
        local_264 = 0xffffffff;
        local_218 = "(i5.5)";
        local_210 = 6;
        local_268 = 0x5000;
        _gfortran_st_write(&local_268);
        _gfortran_transfer_integer_write(&local_268,param_1,4);
        _gfortran_st_write_done(&local_268);
        get_isotropy_filename_
                  (&local_48,param_3,param_4,&local_40,&local_58,param_7,param_8,param_9,param_10);
        if (local_58 == 0) {
          if (local_40 == 0) {
            get_isotropy_(&local_44,param_1,param_3,&DAT_005a602c,param_8,param_9);
          }
        }
        else {
          local_260 = "../source/get_isotropy_file.f";
          local_258 = 0xaf;
          local_218 = "(1x,a)";
          local_210 = 6;
          local_268 = 0x1000;
          local_264 = 6;
          _gfortran_st_write(&local_268);
          _gfortran_transfer_character_write
                    (&local_268,"You have requested information about an isotropy subgroup for: ",
                     0x3f);
          _gfortran_st_write_done(&local_268);
          local_260 = "../source/get_isotropy_file.f";
          local_258 = 0xb2;
          local_218 = "(6x,a7,a8,a14,a8)";
          local_210 = 0x11;
          local_268 = 0x1000;
          local_264 = 6;
          _gfortran_st_write(&local_268);
          _gfortran_transfer_character_write(&local_268,"irrep: ",7);
          _gfortran_transfer_character_write(&local_268,lVar4 * 8 + 0x1509ae6,8);
          _gfortran_transfer_character_write(&local_268," space group: ",0xe);
          _gfortran_transfer_character_write(&local_268,(long)local_44 * 8 + 0xe976fe,8);
          _gfortran_st_write_done(&local_268);
          local_260 = "../source/get_isotropy_file.f";
          local_258 = 0xb4;
          local_218 = "(1x,a)";
          local_210 = 6;
          local_268 = 0x1000;
          local_264 = 6;
          _gfortran_st_write(&local_268);
          _gfortran_transfer_character_write
                    (&local_268,"This data base is already being generated.  Try again later.",0x3c)
          ;
          _gfortran_st_write_done(&local_268);
          *param_6 = 1;
        }
      }
      else {
        if (__isotropy_module_MOD_readfromcmd == 0) {
          local_260 = "../source/get_isotropy_file.f";
          local_258 = 0x9f;
          local_218 = "(1x,a)";
          local_210 = 6;
          local_268 = 0x1000;
          local_264 = 6;
          _gfortran_st_write(&local_268);
          _gfortran_transfer_character_write
                    (&local_268,"Error:  dimension of order parameter for irrep ",0x2f);
          _gfortran_transfer_character_write(&local_268,lVar4 * 8 + 0x1509ae6,8);
          _gfortran_st_write_done(&local_268);
          local_260 = "../source/get_isotropy_file.f";
          local_258 = 0xa1;
          local_218 = "(1x,a)";
          local_210 = 6;
          local_268 = 0x1000;
          local_264 = 6;
          _gfortran_st_write(&local_268);
          _gfortran_transfer_character_write(&local_268,"does not match dimension of irrep",0x21);
          _gfortran_st_write_done(&local_268);
        }
        *param_6 = 2;
      }
    }
    else {
      bomb2_("get_isotropy_file 188",0x15);
    }
  }
  else if (0 < (long)param_10) {
    memset(param_4,0x20,param_10);
  }
  return;
}


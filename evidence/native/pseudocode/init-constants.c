/* DERIVED GHIDRA PSEUDOCODE; NOT ORIGINAL SOURCE; NOT COMPILABLE CODE.
 * Evidence: ev_c3934e6205016da351a9a92ccc97319261a336f2661c83e95d6a8d0a5440cfe6
 * Raw response: evidence/native/rea-linux-iso-recovered/07-init-constants.json
 * Types, array bounds and control flow may be imperfectly recovered.
 */

void init_constants_(void)

{
  size_t sVar1;
  int iVar2;
  void *__ptr;
  ulong *puVar3;
  long lVar4;
  size_t __size;
  uint local_2c8;
  undefined4 local_2c4;
  char *local_2c0;
  undefined4 local_2b8;
  int *local_2a0;
  size_t local_290;
  void *local_288;
  undefined *local_280;
  undefined8 local_278;
  undefined8 local_230;
  char *local_228;
  undefined4 local_198;
  int local_b4;
  undefined4 local_b0;
  int local_ac;
  undefined1 local_a8 [120];
  
  _gfortran_getenv("ISODATA",local_a8,7,100);
  find_io_unit_(&local_b0);
  local_2c0 = "../source/constants.f";
  local_2b8 = 0x13;
  local_ac = 0;
  local_2a0 = &local_ac;
  iVar2 = _gfortran_string_len_trim(100,local_a8);
  lVar4 = (long)iVar2;
  if (lVar4 < 0) {
    lVar4 = 0;
  }
  sVar1 = lVar4 + 9;
  __size = 1;
  if (sVar1 != 0) {
    __size = sVar1;
  }
  __ptr = malloc(__size);
  _gfortran_concat_string(sVar1,__ptr,lVar4,local_a8,9,"const.dat");
  local_280 = &DAT_005a50a3;
  local_278 = 3;
  local_228 = "read";
  local_230 = 4;
  local_198 = 0;
  local_2c8 = 0x1004320;
  local_2c4 = local_b0;
  local_290 = sVar1;
  local_288 = __ptr;
  _gfortran_st_open(&local_2c8);
  free(__ptr);
  if (local_ac != 0) {
    local_2c0 = "../source/constants.f";
    local_2b8 = 0x16;
    local_2c8 = 0x80;
    local_2c4 = 6;
    _gfortran_st_write(&local_2c8);
    _gfortran_transfer_character_write(&local_2c8,"Error: cannot open file ",0x18);
    iVar2 = _gfortran_string_len_trim(100,local_a8);
    lVar4 = (long)iVar2;
    if (lVar4 < 0) {
      lVar4 = 0;
    }
    _gfortran_transfer_character_write(&local_2c8,local_a8,lVar4);
    _gfortran_transfer_character_write(&local_2c8,"const.dat",9);
    _gfortran_st_write_done(&local_2c8);
    bomb2_("constants 23",0xc);
  }
  do {
    local_2c0 = "../source/constants.f";
    local_2b8 = 0x19;
    local_2c8 = 0x88;
    local_2c4 = local_b0;
    _gfortran_st_read(&local_2c8);
    _gfortran_transfer_integer(&local_2c8,&local_b4,4);
    iVar2 = local_b4;
    _gfortran_transfer_real(&local_2c8,&const_common_ + (long)local_b4 * 2,8);
    _gfortran_st_read_done(&local_2c8);
  } while ((local_2c8 & 3) != 2);
  local_2c0 = "../source/constants.f";
  local_2b8 = 0x1b;
  local_2c8 = 0;
  local_2c4 = local_b0;
  _gfortran_st_close(&local_2c8);
  DAT_01bb1244 = iVar2 + 1;
  if (2 < DAT_01bb1244) {
    puVar3 = &DAT_01bb1250;
    do {
      puVar3[1] = *puVar3 ^ 0x8000000000000000;
      puVar3 = puVar3 + 2;
    } while (puVar3 != &DAT_01bb1260 + (ulong)(iVar2 - 2U >> 1) * 2);
  }
  return;
}


#define KG2_DER_POINT \
      dxgxx, dxkxx, cdxkxx, \
      dxgxy, dxkxy, cdxkxy, \
      dxgxz, dxkxz, cdxkxz, \
      dxgyy, dxkyy, cdxkyy, \
      dxgyz, dxkyz, cdxkyz, \
      dxgzz, dxkzz, cdxkzz, \
      dygxx, dykxx, cdykxx, \
      dygxy, dykxy, cdykxy, \
      dygxz, dykxz, cdykxz, \
      dygyy, dykyy, cdykyy, \
      dygyz, dykyz, cdykyz, \
      dygzz, dykzz, cdykzz, \
      dzgxx, dzkxx, cdzkxx, \
      dzgxy, dzkxy, cdzkxy, \
      dzgxz, dzkxz, cdzkxz, \
      dzgyy, dzkyy, cdzkyy, \
      dzgyz, dzkyz, cdzkyz, \
      dzgzz, dzkzz, cdzkzz, \
      d2xx_gxx,\
      d2xx_gxy,\
      d2xx_gxz,\
      d2xx_gyy,\
      d2xx_gyz,\
      d2xx_gzz,\
      d2xy_gxx,\
      d2xy_gxy,\
      d2xy_gxz,\
      d2xy_gyy,\
      d2xy_gyz,\
      d2xy_gzz,\
      d2xz_gxx,\
      d2xz_gxy,\
      d2xz_gxz,\
      d2xz_gyy,\
      d2xz_gyz,\
      d2xz_gzz,\
      d2yy_gxx,\
      d2yy_gxy,\
      d2yy_gxz,\
      d2yy_gyy,\
      d2yy_gyz,\
      d2yy_gzz,\
      d2yz_gxx,\
      d2yz_gxy,\
      d2yz_gxz,\
      d2yz_gyy,\
      d2yz_gyz,\
      d2yz_gzz,\
      d2zz_gxx,\
      d2zz_gxy,\
      d2zz_gxz,\
      d2zz_gyy,\
      d2zz_gyz,\
      d2zz_gzz


#define DECLARE_KG2_DER_POINT \
        CCTK_REAL dxkxx, dxgxx, cdxkxx &&\
        CCTK_REAL dxkxy, dxgxy, cdxkxy &&\
        CCTK_REAL dxkxz, dxgxz, cdxkxz &&\
        CCTK_REAL dxkyy, dxgyy, cdxkyy &&\
        CCTK_REAL dxkyz, dxgyz, cdxkyz &&\
        CCTK_REAL dxkzz, dxgzz, cdxkzz &&\
        CCTK_REAL dykxx, dygxx, cdykxx &&\
        CCTK_REAL dykxy, dygxy, cdykxy &&\
        CCTK_REAL dykxz, dygxz, cdykxz &&\
        CCTK_REAL dykyy, dygyy, cdykyy &&\
        CCTK_REAL dykyz, dygyz, cdykyz &&\
        CCTK_REAL dykzz, dygzz, cdykzz &&\
        CCTK_REAL dzkxx, dzgxx, cdzkxx &&\
        CCTK_REAL dzkxy, dzgxy, cdzkxy &&\
        CCTK_REAL dzkxz, dzgxz, cdzkxz &&\
        CCTK_REAL dzkyy, dzgyy, cdzkyy &&\
        CCTK_REAL dzkyz, dzgyz, cdzkyz &&\
        CCTK_REAL dzkzz, dzgzz, cdzkzz &&\
        CCTK_REAL d2xx_gxx &&\
        CCTK_REAL d2xx_gxy &&\
        CCTK_REAL d2xx_gxz &&\
        CCTK_REAL d2xx_gyy &&\
        CCTK_REAL d2xx_gyz &&\
        CCTK_REAL d2xx_gzz &&\
        CCTK_REAL d2xy_gxx &&\
        CCTK_REAL d2xy_gxy &&\
        CCTK_REAL d2xy_gxz &&\
        CCTK_REAL d2xy_gyy &&\
        CCTK_REAL d2xy_gyz &&\
        CCTK_REAL d2xy_gzz &&\
        CCTK_REAL d2xz_gxx &&\
        CCTK_REAL d2xz_gxy &&\
        CCTK_REAL d2xz_gxz &&\
        CCTK_REAL d2xz_gyy &&\
        CCTK_REAL d2xz_gyz &&\
        CCTK_REAL d2xz_gzz &&\
        CCTK_REAL d2yy_gxx &&\
        CCTK_REAL d2yy_gxy &&\
        CCTK_REAL d2yy_gxz &&\
        CCTK_REAL d2yy_gyy &&\
        CCTK_REAL d2yy_gyz &&\
        CCTK_REAL d2yy_gzz &&\
        CCTK_REAL d2yz_gxx &&\
        CCTK_REAL d2yz_gxy &&\
        CCTK_REAL d2yz_gxz &&\
        CCTK_REAL d2yz_gyy &&\
        CCTK_REAL d2yz_gyz &&\
        CCTK_REAL d2yz_gzz &&\
        CCTK_REAL d2zz_gxx &&\
        CCTK_REAL d2zz_gxy &&\
        CCTK_REAL d2zz_gxz &&\
        CCTK_REAL d2zz_gyy &&\
        CCTK_REAL d2zz_gyz &&\
        CCTK_REAL d2zz_gzz &&\


#include <oxstd.h>


rbinom(const a,const b, const n, const p)

 {
  decl u,sm;

  u = ranu(a,b);

  sm=(u-p)./fabs(u-p);
  sm=(1-sm)/2;
  
  return sm;
 }

#include <oxstd.h>


rtnorm(const r,const c, const mu, const sigma2, const a, const b)

 {
  decl u;

  u = sqrt(sigma2).*quann( probn((a-mu)./sqrt(sigma2))+ pow(10,-100) + ranu(r,c) .* ( probn((b-mu)./sqrt(sigma2)) - probn((a-mu)./sqrt(sigma2)))) + mu;
  return u;
 }

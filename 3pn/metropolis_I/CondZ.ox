CondZ(const a, const b, const c, const theta,const V, const Resp)
{
decl NumInd, NumItens, Z;

Z= Resp .==0 .? 0 .: (V.==1 .? 1 .: 0);

return Z;
}
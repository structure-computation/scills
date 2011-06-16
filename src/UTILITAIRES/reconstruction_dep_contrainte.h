
struct assign_dep_contrainte{
template<class SST> void operator()(SST &S) const{
S.f->get_result()=S.q;
S.f->update_variables();
S.f->call_after_solve();
//S.f.vectors[0]=S.q;
//S.f.update_variables();
//S.f.call_after_solve();
}
};


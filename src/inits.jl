function initialize_steady_vars(geo::LHPGeometry)
    return ComponentArray(
        grv     = (T=zeros(geo.grv.N),  P=zeros(geo.grv.N),     x=ones(geo.grv.N),      mass=zeros(geo.grv.N)),
        vl      = (T=zeros(geo.vl.N),   P=zeros(geo.vl.N),      x=ones(geo.vl.N),       mass=zeros(geo.vl.N)),
        cond    = (T=zeros(geo.cond.N), P=zeros(geo.cond.N),    x=ones(geo.cond.N),     mass=zeros(geo.cond.N)),
        ll      = (T=zeros(geo.ll.N),   P=zeros(geo.ll.N),      x=zeros(geo.ll.N),      mass=zeros(geo.ll.N)),
        bay     = (T=zeros(geo.bay.N),  P=zeros(geo.bay.N),     x=zeros(geo.bay.N),     mass=zeros(geo.bay.N)),
        cc      = (T=zeros(geo.cc.N),   P=zeros(geo.cc.N),      x=zeros(geo.cc.N),      mass=zeros(geo.cc.N))
    )
end
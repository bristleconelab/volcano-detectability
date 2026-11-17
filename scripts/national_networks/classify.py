import json,csv,math,re,sys,collections
D=sys.argv[1]; OUT=sys.argv[2]
def hav(a,b,c,d):
    a,b,c,d=map(math.radians,(a,b,c,d)); h=math.sin((c-a)/2)**2+math.cos(a)*math.cos(c)*math.sin((d-b)/2)**2
    return 6371*2*math.asin(math.sqrt(h))
OPEN={'NGL','UNAVCO','PBO','USGS'}
rows=list(csv.DictReader(open(D+'/the70.csv')))
world=json.load(open(D+'/dl_def.json'))['deformation']
# instrument types by (code) gathered from all 70 per-volcano files
itype=collections.defaultdict(set)
for r in rows:
    j=json.load(open(D+'/vol/%s.json'%r['vnum'])); de=j['deformation']
    code={s['station_id']:s['station_code'].strip() for s in de['station']}
    for g in ('general','tilt'):
        for i in de['instrument'][g]:
            if i['station_id'] in code: itype[code[i['station_id']]].add(i['inst_type_abbre'])
def key(code):
    c=code.strip()
    m=re.match(r'NGL_(.{4})_\d+$',c) or re.match(r'UNAVCO_(.{4})$',c)
    if m: return m.group(1).rstrip('_')
    m=re.match(r'(.+)_GPS$',c)
    if m: return m.group(1)
    return c
def isgnss(code,ty):
    return bool(re.match(r'(NGL|UNAVCO)_',code)) or code.endswith('_GPS') or any('GPS' in t or 'GNSS' in t for t in ty)
out=[]
for r in rows:
    la,lo=float(r['lat']),float(r['lon'])
    j=json.load(open(D+'/vol/%s.json'%r['vnum']))
    cand=[]
    for s in j['deformation']['station']:
        src={s.get('datasource_%d'%k) for k in (1,2,3)}-{None}
        cand.append((s['station_code'].strip(),float(s['station_lat']),float(s['station_lon']),src))
    for w in world:
        src={w['datasource1'],w['datasource2'],w['datasource3']}-{None}
        cand.append((w['station_code'].strip(),float(w['station_lat']),float(w['station_lon']),src))
    g_open,g_nat,oth=collections.defaultdict(set),collections.defaultdict(set),collections.defaultdict(set)
    for code,a,b,src in cand:
        if hav(la,lo,a,b)>30: continue
        ty=itype.get(code,set()); k=key(code)
        if isgnss(code,ty):
            (g_open if src<=OPEN or src&{'NGL','UNAVCO'} and not src-OPEN-{'SNET'} else g_nat)[k]|=src
        else:
            oth[k]|=src|{'type:'+'/'.join(sorted(ty)) if ty else 'type:unrecorded'}
    out.append(dict(vnum=r['vnum'],name=r['name'],country=r['country'],
        gnss_open=sorted(g_open),gnss_nat=sorted(g_nat),
        nat_src=sorted(set().union(*g_nat.values())) if g_nat else [],
        other=len(oth),other_types=sorted(set(x for v in oth.values() for x in v if x.startswith('type:'))),
        other_src=sorted(set(x for v in oth.values() for x in v if not x.startswith('type:'))),
        seis=len(j['seismic']['station'])))
json.dump(out,open(OUT,'w'),indent=0)
for o in out: print(o['vnum'],o['name'][:16].ljust(16),'open',o['gnss_open'],'nat',o['gnss_nat'],o['nat_src'],'oth',o['other'],o['other_types'],o['other_src'])

import json,csv,sys
D=sys.argv[1]; OUT=sys.argv[2]
cl={c['vnum']:c for c in json.load(open(D+'/classified.json'))}
G='https://wovodat.org/gvmid/index.php?type=single&vd_num='
JDR='https://www.jstage.jst.go.jp/article/jdr/14/1/14_6/_article'
# lit: vnum -> (kind, text, public, url, verified)  kind: cont=continuous GNSS, camp=campaign/unspecified GNSS, oth=non-GNSS deformation only
L={
'261080':('cont','CVGHM cGPS: 4 stations since Feb 2011 (Hendrasto 2012 JDR); CVGHM cGPS KBYK/LKWR/MRDG/SNBG used by Siregar 2018; Hotta 2018 JVGR models 2013-2016 cGPS','no','https://core.ac.uk/download/211044799.pdf','search-summary'),
'263280':('cont','SATREPS seismic+GNSS network installed Dec 2014-Sep 2015 at Kelud, data to CVGHM decision system (Nakamichi 2019 JDR); 2 tiltmeters before 2014 eruption (JpGU 2015)','no',JDR,'page-read (JDR); tilt from search-summary'),
'263300':('cont','SATREPS GNSS 2014-2015 (Nakamichi 2019 JDR); PVMBG CORS stations (CLUM, LEKR, SAWR) 2015-2023 used in ITS thesis','no',JDR+' ; https://repository.its.ac.id/109307','page-read'),
'263140':('cont','SATREPS seismic+GNSS network at Galunggung, Dec 2014-Sep 2015 (Nakamichi 2019 JDR)','no',JDR,'page-read'),
'264020':('cont','continuous GPS around Agung used for 2018-2019 eruption (AGU 2020 abstract V016-0008)','no','https://ui.adsabs.harvard.edu/abs/2020AGUFMV016.0008G/abstract','search-summary'),
'266100':('camp','GPS observation points at Lokon 2009-2011, 1.1-7 cm/yr (IJOG; CVGHM+ITB authors; continuous vs campaign not stated)','no','https://ijog.geologi.esdm.go.id/index.php/IJOG/article/view/147','page-read'),
'264010':('camp','campaign GPS network of 11 points (Kimata 2001 report); periodic GPS 2008/2009/2013/2015 (Faris 2018)','no','https://www2.jpgu.org/meeting/2001/pdf/jp/jp-007_e.pdf','search-summary'),
'263350':('camp','GPS surveys 2002-2005 (dm-level shifts); tiltmeter/EDM monitoring described','no','https://core.ac.uk/download/426986413.pdf','search-summary'),
'263100':('camp','GPS surveys before/during/after 2002 eruption (Abidin et al., ITB)','no','https://lppm.itb.ac.id/?p=49715','search-summary'),
'273010':('cont','PHIVOLCS bulletin: real-time electronic tilt + continuous GPS monitoring','no','https://www.phivolcs.dost.gov.ph/bulusan-volcano-bulletin-19-may-2021-08-00-a-m/','page-read'),
'341040':('cont','UNAM permanent 5-station network with Trimble NetR9 GPS at 2.2-2.4 km from crater (EGU2011-4088, Cabral-Cano et al.)','no','https://meetingorganizer.copernicus.org/EGU2011/EGU2011-4088.pdf','page-read'),
'351020':('cont','SGC Sep 2026 bulletin: deformation monitored with electronic tiltmeters, GNSS stations and DInSAR','no','https://www2.sgc.gov.co/Noticias/boletinesDocumentos/Boletin_Informativo_No_3020_septiembre_2026.pdf','page-read'),
'351060':('cont','SGC permanent GNSS since 2009; 20 GNSS stations across Ruiz, Machin, Purace, Galeras (2015 paper)','no','https://iris.uniroma1.it/handle/11573/779663','search-summary'),
'351080':('cont','SGC permanent GNSS since 2009; 20 GNSS stations across Ruiz, Machin, Purace, Galeras (2015 paper)','no','https://iris.uniroma1.it/handle/11573/779663','search-summary'),
'345070':('cont','OVSICORI-UNA: 2 CGPS at Turrialba since 2010 (COCONet slides); Turrialba-Irazu cGPS network 2014-2016 (2024 study)','no','https://coconet.unavco.org/lib/downloads/COCONet_Costa%20Rica_protti.pdf','search-summary'),
'345060':('cont','OVSICORI-UNA: semipermanent GPS at Irazu since 1999 (COCONet slides); Turrialba-Irazu cGPS network 2014-2016','no','https://coconet.unavco.org/lib/downloads/COCONet_Costa%20Rica_protti.pdf','search-summary'),
'345040':('oth','Poas deformation network: precise levelling, geodetic and EDM networks (UNA technical paper); OVSICORI cGPS at Poas not verified','unknown','https://repositorio.una.ac.cr/handle/11056/23913','search-summary'),
'343100':('cont','cGPS VMIG installed 2007 by MARN + Univ. Wisconsin (archived at UNAVCO, now inactive in NGL); MARN-INGV GPS after 2013 eruption','yes','https://oldri.ues.edu.sv/id/eprint/12791/1/19201048.pdf','search-summary'),
'343020':('cont','cGPS SNJE installed 2007 by MARN + Univ. Wisconsin (UNAVCO archive, inactive in NGL); 13-station campaign pilot (MTU)','yes','https://digitalcommons.mtu.edu/geo-fp/4','search-summary'),
'343030':('cont','cGPS SNJE (7.2 km) installed 2007 by MARN + Univ. Wisconsin (UNAVCO archive, inactive in NGL); same station as Santa Ana','yes','https://digitalcommons.mtu.edu/geo-fp/4','search-summary'),
'360150':('cont','UWI-SRC: cGPS SVGB + 4 cGPS added (SVGR, SVGS, SVGF, SVGG) by Feb 2021 (Joseph et al. 2022 Nat Commun)','unknown','https://www.nature.com/articles/s41467-022-31901-4','search-summary'),
'252140':('cont','RVO/Geoscience Australia continuous GNSS at Rabaul caldera since ~2000 (Saunders et al. 2023 Geosciences)','no','https://www.mdpi.com/2076-3263/13/8/249','search-summary'),
'275001':('camp','Hainan Earthquake Agency volcano network since 2007; 7 GNSS networks incl. multi-period mobile GNSS (campaign) in Qiongbei area','no','https://www.dzdz.ac.cn/CN/Y2022/V44/I5/1128','page-read'),
'342030':('oth','INSIVUMEH Santiaguito: electronic tiltmeter, dry tilt, levelling lines (WOVO directory); no GNSS found','unknown','https://wovo.iavceivolcano.org/members-observatories-directory/guatemala-observatorio-vulcanologica-del-volcan-santiaguito-insivumeh/','search-summary'),
}
rows=[]; summ={'any_gnss':[], 'cont_nat':[], 'only_hist_open':[], 'camp_only':[], 'oth_only':[], 'none':[]}
for v,c in cl.items():
    gn=sorted(set(c['gnss_open'])|set(c['gnss_nat']))
    natsrc=[s for s in c['nat_src'] if s not in ('NGL','UNAVCO','PBO','USGS')]
    other=c['other'] - len(set(gn) & set())  # other already excludes typed GNSS
    ev=[]
    if gn:
        # station codes and locations from GVMID are not redistributed (WOVOdat data policy); counts only
        ev.append('GVMID<=30km GNSS stations: %d'%len(gn))
        if c['gnss_open']: ev.append('of which open-archive (NGL/UNAVCO) sites: %d, none active in current NGL list'%len(c['gnss_open']))
        if c['gnss_nat']: ev.append('of which observatory GNSS (%s): %d'%('/'.join(natsrc),len(c['gnss_nat'])))
    if c['other']: ev.append('GVMID other deformation %d (%s; %s)'%(c['other'],'/'.join(t[5:] for t in c['other_types']),'/'.join(c['other_src'])))
    if not gn and not c['other']: ev.append('GVMID: no deformation station <=30 km')
    lit=L.get(v)
    if lit: ev.append('LIT [%s]: %s'%(lit[4],lit[1]))
    # public
    if c['gnss_nat'] or (lit and lit[2]=='no'): pub='no'
    elif lit and lit[2]=='unknown': pub='unknown'
    elif gn: pub='yes'
    elif c['other']: pub='no'
    else: pub='unknown'
    url=G+v + (' ; '+lit[3] if lit else '')
    # categories
    has_cont_nat = bool(c['gnss_nat']) or (lit and lit[0]=='cont')
    if gn or (lit and lit[0] in ('cont','camp')): summ['any_gnss'].append(c['name'])
    if has_cont_nat: summ['cont_nat'].append(c['name'])
    elif gn and not (lit and lit[0]=='camp'): summ['only_hist_open'].append(c['name'])
    elif lit and lit[0]=='camp': summ['camp_only'].append(c['name'])
    elif c['other'] or (lit and lit[0]=='oth'): summ['oth_only'].append(c['name'])
    else: summ['none'].append(c['name'])
    cat=('observatory_cgnss_not_open' if has_cont_nat else 'open_sites_inactive' if gn and not (lit and lit[0]=='camp')
         else 'campaign_gps_only' if lit and lit[0]=='camp' else 'other_deformation_only' if c['other'] or (lit and lit[0]=='oth')
         else 'none_found')
    rows.append(dict(vnum=v,name=c['name'],country=c['country'],category=cat,gvmid_gnss_n=len(gn),gvmid_other_deformation_n=c['other'],
        evidence='; '.join(ev),data_public=pub,source_url=url))
with open(OUT,'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
summ['obs_nonGNSS_deformation']=[c['name'] for c in cl.values() if c['other'] and set(c['other_src'])-{'UNAVCO','NGL','PBO','USGS','SNET'}]
for k,v in summ.items(): print(k,len(v),v)

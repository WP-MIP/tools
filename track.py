#!/usr/bin/env python3
import pandas as pd
import rpnpy.librmn.all as rmn
from datetime import datetime as dt
import numpy as np
import os, math

class storm(dict):
    
    def __init__(self, df, id, dpath):
        self.id = id
        self.track = df.loc[self.id]
        self.num = int(self.track['NUMBER'].iloc[0])
        self.dpath = dpath
        self.start_str = self.track['ISO_TIME'].iloc[0]
        self.start = dt.fromisoformat(self.start_str)
        self.end_str = self.track['ISO_TIME'].iloc[-1]
        self.end = dt.fromisoformat(self.end_str)
        self.dxy = 8
        self._find_ics()
    def tctrack(self):
        fd = rmn.fstopenall(os.path.join(self.dpath, 'gh_1000*'))
        fcst_keys = {}
        # Retrieve list of keys for each initialization
        for ic in self.ic_dates:
            fcst_keys[ic] = []
            keys = rmn.fstinl(fd)
            for key in keys:
                meta = rmn.fstprm(key)
                if meta['dateo'] == ic:
                    if rmn.difdatr(meta['datev'], self._torpndate(self.start)) >= 0 \
                       and rmn.difdatr(meta['datev'], self._torpndate(self.end)) <= 0:
                        fcst_keys[ic].append(key)
        # Process each initialization
        self.tracks = {}
        for ic in fcst_keys.keys():
            # Identify initial storm location
            self.tracks[ic] = []
            for key in fcst_keys[ic]:
                fld = rmn.fstluk(key)
                grid = rmn.readGrid(fd, fld)
                (lat, lon, date, fcst) = self._btPos(key)
                found = self._findCentre(grid, fld, lat, lon)
                found['date'] = date
                found['fcst'] = fcst
                (found['pres'], found['wind']) = self._vitals(found)
                self.tracks[ic].append(found)
        rmn.fstcloseall(fd)
        return()
    def write(self, path):
        for ic in self.tracks.keys():
            idat = self._todtime(ic)
            fname = os.path.join(path, '0001'+idat.strftime('%Y%m%d%H')+'_FC_000360_enp')
            try:
                os.makedirs(path)
            except FileExistsError:
                pass
            fd = open(fname, "a")
            fd.write("00000 {0} M={1:2d}   1 SNBR={2:4d}\n".format(idat.strftime('%d/%m/%Y'), len(self.tracks[ic]),
                                                                   self.num))
            ext = '*0000000*00000000000000000000*00000000000000000000*00000000000000000000*\n'
            for fix in self.tracks[ic]:
                wind = fix['min'] > 0 and round(fix['wind']) or 0
                fd.write("00000 {0}*{1}{2:3d}{3:4d}{4:5d}".format(fix['date'].strftime('%Y/%m/%d/%H'),round(fix['lat']*10),
                                                                  round(fix['lon']*10), wind, round(fix['pres'])) + ext)
            fd.write("00000  EX\n")
            fd.close()
    def _find_ics(self):
        fd = rmn.fstopenall(os.path.join(self.dpath, 'gh_1000*'))
        end = self._torpndate(self.end)
        keys = rmn.fstinl(fd, datev=end)
        self.ic_dates = []
        for key in keys:
            meta = rmn.fstprm(key)
            self.ic_dates.append(meta['dateo'])
        rmn.fstcloseall(fd)
        return()
    def _btPos(self, key):
        meta = rmn.fstprm(key)
        irow = self.track.loc[self.track['ISO_TIME'] == self._rpntoibtracs(meta['datev'])]
        (ilat, ilon) = (float(irow['LAT'].iloc[0]), float(irow['LON'].iloc[0]))
        if ilon < 0: ilon = ilon + 360
        return((ilat, ilon, self._todtime(meta['datev']), meta['ip2']))
    def _findCentre(self, grid, fld, lat, lon):
        xy = rmn.gdxyfll(grid['id'], lat, lon)
        x = int(xy['x'][0])-1
        y = int(xy['y'][0])-1
        search_fld = fld['d'][x-self.dxy:x+self.dxy+1,y-self.dxy:y+self.dxy+1]
        xy_sub = np.where(search_fld == search_fld.min())
        ll = rmn.gdllfxy(grid['id'], [x-self.dxy+xy_sub[0]], [y-self.dxy+xy_sub[1]])
        xf = int(ll['x'][0,0])
        yf = int(ll['y'][0,0])
        fldf = float(fld['d'][xf,yf])
        ismin = 0
        if fld['d'][max(xf-1,0):xf+2,yf-1:yf+2].min() == fldf:
            ismin=1
        return({'x':xf, 'y':yf, 'hght':fldf, 'min':ismin, 
                'lat':float(ll['lat'][0,0]), 'lon':float(ll['lon'][0,0])})
    def _vitals(self, fix):
        wind_dxy = 4
        x = fix['x']
        y = fix['y']
        fldt = self._readfld(fix, 't')
        pres = 1000 * math.exp(9.81 * fix['hght'] / (287. * float(fldt['d'][x,y])))
        fldu = self._readfld(fix, 'u')
        fldv = self._readfld(fix, 'v')
        wspd = np.sqrt(fldu['d']**2 + fldv['d']**2)
        wind = wspd[max(x-wind_dxy,0):x+wind_dxy+1,y-wind_dxy:y+wind_dxy+1].max()
        return(pres, wind)
    def _readfld(self, fix, name):
        fdfld = rmn.fstopenall(os.path.join(self.dpath, name+'_1000*'))
        fld = rmn.fstlir(fdfld, datev=self._torpndate(fix['date']), ip2=fix['fcst'])
        rmn.fstcloseall(fdfld)
        return(fld)
    def _torpndate(self, dtime):
        return(rmn.newdate(rmn.NEWDATE_PRINT2STAMP, int(dtime.strftime('%Y%m%d')),
                           int(dtime.strftime('%H%M%S'))*100))
    def _todtime(self, rpndate):
        (ymd, hms) = rmn.newdate(rmn.NEWDATE_STAMP2PRINT, rpndate)
        return(dt.strptime(str(ymd)+str(int(hms/1000000)).rjust(2, '0'), '%Y%m%d%H'))
    def _dttoibtracs(self, dtime):
        return(dtime.strftime('%Y-%m-%d %H:%M:%S'))
    def _rpntoibtracs(self, rpndate):
        return(self._dttoibtracs(self._todtime(rpndate)))
        

if __name__ == "__main__":

    tc_wp = ["2024141N03142", "2024151N18113", "2024224N27154", "2024225N22135",
             "2024225N24147", "2024231N24126", "2024246N22147", "2024259N12145",
             "2024267N29129", "2024269N14150", "2024278N11150", "2024298N13150"]
    tc_na = ["2024181N09320", "2024216N20284", "2024225N14313", "2024253N21266",
             "2024269N39302", "2024274N14328", "2024279N21265"]
    
    fcst_path = "data/oic/cwao/pm"
    opath = "tracks/oic/cwao/pm"

    # North Pacific tropical cyclone tracking
    df_wp = pd.read_csv('ibtracs/ibtracs.WP.list.v04r01.csv', index_col="SID")
    for tc in tc_wp:
        bt = storm(df_wp, tc, fcst_path)
        bt.tctrack()
        bt.write(os.path.join(opath,'WP'))

    # North Atlantic tropical cyclone tracking
    df_na = pd.read_csv('ibtracs/ibtracs.NA.list.v04r01.csv', index_col="SID")
    for tc in tc_na:
        bt = storm(df_na, tc, fcst_path)
        bt.tctrack()
        bt.write(os.path.join(opath,'NA'))
    


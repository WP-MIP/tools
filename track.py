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
        self.search_rad = 8
        self._find_ics()
    def tctrack(self):
        fd = rmn.fstopenall(os.path.join(self.dpath, 'gh_1000*'))
        fcst_keys = {}
        et = {}
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
                        nature = self.track.loc[self.track['ISO_TIME'] ==
                                                self._rpntoibtracs(meta['datev'])]['NATURE'].iloc[0]
                        if ic not in et.keys() and (nature == "ET" or key == keys[-1]):
                                et[ic] = {'date':self._todtime(meta['datev']), 'key':key}
        # Process each initialization
        self.tracks = {}
        for ic in fcst_keys.keys():
            # Identify initial storm location
            self.tracks[ic] = []
            ikey0 = fcst_keys[ic].index(et[ic]['key'])
            track_fwd = self._tracking(fd, fcst_keys[ic][ikey0:len(fcst_keys[ic])+1], 1)
            track_back = self._tracking(fd, fcst_keys[ic][0:ikey0+1], -1)
            track = track_back + track_fwd[1:]
            if len(track) > 1:
                self.tracks[ic].append(track_back + track_fwd[1:])
        rmn.fstcloseall(fd)
        return()
    def write(self, path):
        for ic in self.tracks.keys():
            if (len(self.tracks[ic]) < 1): continue
            idat = self._todtime(ic)
            fname = os.path.join(path, '0001'+idat.strftime('%Y%m%d%H')+'_FC_000360_enp')
            try:
                os.makedirs(path)
            except FileExistsError:
                pass
            try:
                with open(fname, "r") as fd:
                    cnt = len(fd.readlines())
            except FileNotFoundError:
                cnt = 0
            with open(fname, "a+") as fd:
                fd.write("{0:05d} {1} M={2:2d}   1 SNBR={3:4d}\n".format(cnt*10, idat.strftime('%d/%m/%Y'),
                                                                         len(self.tracks[ic]), self.num))
                ext = '*0000000*00000000000000000000*00000000000000000000*00000000000000000000*\n'
                for fixes in self.tracks[ic]:
                    for fix in fixes:
                        wind = fix['min'] > 0 and round(fix['wind']) or 0
                        cnt = cnt+1
                        fixstr = "{0:05d} {1}*{2:03d}{3:04d}{4:4d}{5:5d}*{6:03d}{7:04d}"
                        fd.write(fixstr.format(cnt*10, fix['date'].strftime('%Y/%m/%d/%H'),
                                               round(fix['lat']*10), round(fix['lon']*10),
                                               wind, round(fix['pres']),
                                               round(fix['latwind']*10), round(fix['lonwind']*10)) + ext)
                fd.write("{0:05d}  EX\n".format(cnt*10+5))
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
    def _fcstPos(self, key, lats, lons, motion=True):
        meta = rmn.fstprm(key)
        datev = self._todtime(meta['datev'])
        ip2 = meta['ip2']
        if motion:
            return((2*lats[-1]-lats[-2], 2*lons[-1]-lons[-2], datev, ip2))
        else:
            return((lats[-1], lons[-1], datev, ip2))
        
    def _findCentre(self, grid, fld, lat, lon, radmul=1):
        xy = rmn.gdxyfll(grid['id'], lat, lon)
        x = int(xy['x'][0])-1
        y = int(xy['y'][0])-1
        rad = round(radmul * self.search_rad)
        search_fld = fld['d'][max(x-rad,0):min(x+rad+1,grid['ni']-1),y-rad:y+rad+1]
        xy_sub = np.where(search_fld == search_fld.min())
        ll = rmn.gdllfxy(grid['id'], [x-rad+xy_sub[0]], [y-rad+xy_sub[1]])
        xf = int(ll['x'][0,0])
        yf = int(ll['y'][0,0])
        fldf = float(fld['d'][xf,yf])
        ismin = 0
        xflim = min(max(xf, 2), grid['ni']-3)
        if fld['d'][xflim-1:xflim+2,yf-1:yf+2].min() == fldf:
            ismin=1
        return({'x':xf, 'y':yf, 'hght':fldf, 'min':ismin, 
                'lat':float(ll['lat'][0,0]), 'lon':float(ll['lon'][0,0])})
    def  _tracking(self, fd, keys_in, dir):
        keys = keys_in[::dir]        
        lats = []
        lons = []
        track = []
        for key in keys:
            fld = rmn.fstluk(key)
            grid = rmn.readGrid(fd, fld)
            if len(lats) < 2:
                (lat, lon, date, fcst) = self._btPos(key)
                found = self._findCentre(grid, fld, lat, lon)
                if found['min'] < 1:
                    found = self._findCentre(grid, fld, lat, lon, radmul=2)
            else:
                (lat, lon, date, fcst) = self._fcstPos(key, lats, lons)
                found = self._findCentre(grid, fld, lat, lon)
                if found['min'] < 1:
                    found = self._findCentre(grid, fld, lat, lon, radmul=2)
                if found['min'] < 1:
                    (lat, lon, date, fcst) = self._fcstPos(key, lats, lons, motion=False)
                    found = self._findCentre(grid, fld, lat, lon)
                if found['min'] < 1:
                    found = self._findCentre(grid, fld, lat, lon, radmul=2)
            if found['min'] < 1:
                break
            lats.append(found['lat'])
            lons.append(found['lon'])
            found['date'] = date
            found['fcst'] = fcst
            (found['pres'], found['wind'], found['latwind'], found['lonwind']) = self._vitals(found, grid)
            track.append(found)
        return(track[::dir])
    def _vitals(self, fix, grid):
        wind_dxy = 6
        x = fix['x']
        y = fix['y']
        fldt = self._readfld(fix, 't')
        pres = 1000 * math.exp(9.81 * fix['hght'] / (287. * float(fldt['d'][x,y])))
        fldu = self._readfld(fix, 'u')
        fldv = self._readfld(fix, 'v')
        wspd = np.sqrt(fldu['d']**2 + fldv['d']**2)
        wind_local = wspd[max(x-wind_dxy,0):min(x+wind_dxy+1,grid['ni']-1),y-wind_dxy:y+wind_dxy+1]
        wind = wind_local.max()
        xy_local = np.where(wind_local == wind)
        ll = rmn.gdllfxy(grid['id'], [x-wind_dxy+xy_local[0]], [y-wind_dxy+xy_local[1]])
        return(pres, wind, ll['lat'][0,0], ll['lon'][0,0])
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

    tcid = {'WP':["2024141N03142", "2024151N18113", "2024224N27154", "2024225N22135",
                  "2024225N24147", "2024231N24126", "2024246N22147", "2024259N12145",
                  "2024267N29129", "2024269N14150", "2024278N11150", "2024298N13150"],
            'NA':["2024181N09320", "2024216N20284", "2024225N14313", "2024253N21266",
                  "2024269N39302", "2024274N14328", "2024279N21265"]}


    #tcid = {'NA':["2024181N09320", "2024216N20284", "2024225N14313", "2024253N21266",
    #              "2024269N39302", "2024274N14328", "2024279N21265"]}
    
    fcst_path = "data/oic/cwao/pm"
    opath = "tracks/oic/cwao/pm"

    # Track identified storms in appropriate basins
    for basin in tcid.keys():
        df = pd.read_csv('ibtracs/ibtracs.'+basin+'.list.v04r01.csv', index_col="SID",
                         low_memory=False)
        for tc in tcid[basin]:
            bt = storm(df, tc, fcst_path)
            bt.tctrack()
            bt.write(os.path.join(opath, basin))

    


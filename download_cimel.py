import datetime as dt
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import os

## Where to find documentation about web requests to aeronet
## https://aeronet.gsfc.nasa.gov/print_web_data_help_v3_new.html
## Request ----> 'wget --no-check-certificate  -q  -O test.out "https://aeronet.gsfc.nasa.gov/cgi-bin/print_web_data_v3?site=La_Paz&year=2024&month=6&day=5&year2=2024&month2=6&day2=6&AOD15=1&AVG=10"


### DOWNLOADING LAST N DAYS
days_to_dwnload = 3
today = dt.datetime.now()
date_start = (today + dt.timedelta(days = -days_to_dwnload)).strftime('%Y-%m-%d')
date_end = today.strftime('%Y-%m-%d')


def get_date_info(d_str):
    d = dt.datetime.strptime(d_str, '%Y-%m-%d')
    dstr = d.strftime('%Y%m%d')
    return d.year, d.month, d.day, dstr

y0,m0,d0,d0str = get_date_info(date_start)
yf,mf,df,dfstr = get_date_info(date_end)

instruments = ['Mount_Chacaltaya', 'La_Paz', 'SANTA_CRUZ_UTEPSA']
codenames = ['Chacaltaya', 'La Paz', 'Santa Cruz']
colors = ['g','b','r']

print('Downloading data...')

for instrument in instruments:
    req = f'wget --no-check-certificate  -q  -O cimel_{instrument}.csv "https://aeronet.gsfc.nasa.gov/cgi-bin/print_web_data_v3?site={instrument}&year={y0}&month={m0}&day={d0}&year2={yf}&month2={mf}&day2={df}&AOD15=1&AVG=10&if_no_html=1"'
    os.system(req)
    print(f'cimel_{instrument}.csv DOWNLOADED')

fnames = [f'cimel_{i}.csv' for i in instruments]


def proc_csv(ifile):
    df2 = pd.read_csv(ifile, skiprows = 5)
    datetime2 = pd.to_datetime(df2['Date(dd:mm:yyyy)'] + ' ' + df2['Time(hh:mm:ss)'], format='%d:%m:%Y %H:%M:%S')
    df2['datetime'] = datetime2
    df2 = df2.loc[:,['datetime', 'AOD_500nm', 'AOD_440nm']]
    return df2

dfs = [proc_csv(f'cimel_{i}.csv') for i in instruments]

def plot_aod500(ax, dfs,places, colors):
    for df, p, c in zip(dfs, places, colors):
        ax.scatter(df.datetime, df[f'AOD_500nm'], marker = 'x', label = f'{p} 500nm', c = c)
        #ax.scatter(df.datetime, df.AOD_440nm, marker = '+', label = f'{p} 440nm', c = c)

def plot_aod440(ax, dfs,places, colors):
    for df, p, c in zip(dfs, places, colors):
        ax.scatter(df.datetime, df[f'AOD_440nm'], marker = '+', label = f'{p} 440nm', c = c)
        #ax.scatter(df.datetime, df.AOD_440nm, marker = '+', label = f'{p} 440nm', c = c)

fig, ax = plt.subplots(figsize=(7,3.5), dpi = 200)
plot_aod500(ax, dfs, codenames, colors)
fig.subplots_adjust(right=0.81, bottom = 0.2)
fig.legend(loc = 'center right', fontsize = 6)
ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
ax.xaxis.set_major_locator(mdates.MonthLocator())
ax.xaxis.set_minor_formatter(mdates.DateFormatter('%d'))
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax.xaxis.set_tick_params(which='major', pad=10)
ax.set_xlabel('Date')
ax.grid(alpha = 0.3)
ax.grid(which='minor', axis = 'x', alpha = 0.3)
ax.set_ylabel('Aerosol Optical Depth')
todayformatted = dt.datetime.strftime(today, '%Y-%m-%d %H:%M UTC')
ax.set_title(f'AOD last {days_to_dwnload} days (last update {todayformatted})', fontsize = 10)
fig.savefig('aod500.png', dpi = 120)



fig, ax = plt.subplots(figsize=(7,3.5), dpi = 200)
plot_aod440(ax, dfs, codenames, colors)
fig.subplots_adjust(right=0.81, bottom = 0.2)
fig.legend(loc = 'center right', fontsize = 6)
ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
ax.xaxis.set_major_locator(mdates.MonthLocator())
ax.xaxis.set_minor_formatter(mdates.DateFormatter('%d'))
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax.xaxis.set_tick_params(which='major', pad=10)
ax.set_xlabel('Date')
ax.grid(alpha = 0.3)
ax.grid(which='minor', axis = 'x', alpha = 0.3)
ax.set_ylabel('Aerosol Optical Depth')
todayformatted = dt.datetime.strftime(today, '%Y-%m-%d %H:%M UTC')
ax.set_title(f'AOD last {days_to_dwnload} days (last update {todayformatted})', fontsize = 10)
fig.savefig('aod440.png', dpi = 120)


for f in fnames: os.remove(f)

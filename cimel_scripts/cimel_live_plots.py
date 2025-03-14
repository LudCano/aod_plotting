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
days_to_dwnload = 7
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
codenames2 = ['chc','lpz','scz']
colors = ['turquoise','b','r']
colors2 = ['darkcyan', 'darkblue', 'darkred']

print('Downloading data...')

for instrument in instruments:
    req = f'wget --no-check-certificate  -q  -O cimel_{instrument}.csv "https://aeronet.gsfc.nasa.gov/cgi-bin/print_web_data_v3?site={instrument}&year={y0}&month={m0}&day={d0}&year2={yf}&month2={mf}&day2={df}&AOD15=1&AVG=10&if_no_html=1"'
    os.system(req)
    print(f'cimel_{instrument}.csv DOWNLOADED')

fnames = [f'cimel_{i}.csv' for i in instruments]




def proc_csv2(ifile):
    df2 = pd.read_csv(ifile, skiprows = 5)
    datetime2 = pd.to_datetime(df2['Date(dd:mm:yyyy)'] + ' ' + df2['Time(hh:mm:ss)'], format='%d:%m:%Y %H:%M:%S')
    df2['datetime'] = datetime2
    df2 = df2.loc[:,['datetime', 'AOD_500nm', 'AOD_440nm']]
    df2 = df2[(df2.AOD_500nm > -1) & (df2.AOD_440nm > -1)]
    return df2

def plot_aod500_day(ax, dfs,places, colors, df2, places2, colors2):
    for df, p, c, p2, c2 in zip(dfs, places, colors, places2, colors2):
        ax.scatter(df.datetime, df[f'AOD_500nm'], marker = 'x', label = f'{p} 500nm', c = c)
        ax.scatter(df2.datetime, df2[f'{p2}'], marker = 'o', label = f'{p} 550nm', c = c2, zorder = 0)



goes_aod = pd.read_csv('goes_scripts/aod_places.csv', parse_dates = ['datetime'])


dfs2 = [proc_csv2(f'cimel_{i}.csv') for i in instruments]

os.remove('current_plots/live_aod.png')

fig, ax = plt.subplots(figsize=(7,3.5), dpi = 200)
plot_aod500_day(ax, dfs2, codenames, colors, goes_aod, codenames2, colors2)
x0 = dt.datetime.now().replace(hour = 10, minute=0, second=0, microsecond=0)
xf = x0 + dt.timedelta(hours = 12)
ax.set_xlim(x0, xf)
ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H'))
fig.subplots_adjust(right=0.81, bottom = 0.2)
fig.legend(loc = 'center right', fontsize = 6)
ax.grid(alpha = 0.3)
ax.set_xlabel('Hour (UTC)')
ax.set_ylabel('Aerosol Optical Depth')
todayformatted = dt.datetime.strftime(today, '%Y-%m-%d %H:%M UTC')
ax.set_title(f'AOD LIVE (last update {todayformatted})', fontsize = 10)
print('FIGURE UPDATED')
fig.savefig('current_plots/live_aod.png', dpi = 120)


for f in fnames: os.remove(f)

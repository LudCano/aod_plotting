# Automatic generation of AOD products

**LIVE** plots can be seen at https://ludcano.github.io/aod_plotting/

This code generates automatically the following products: 
   1. AOD from the CIMEL instrument located at three different stations: La Paz, Chacaltaya and Santa Cruz **every 30 minutes** including GOES data for the same timeframe (*)
   2. Composite image for GOES AOD in the South America Region **every day at 15:00 UTC** (*) is displayed in the visualizer and saved in the folder (composites)[https://github.com/LudCano/aod_plotting/tree/main/composites]
   3. Last week AOD (daily mean) for the same CIMEL instruments **every day at 12:00 UTC**
   4. Daily AOD from CIMEL instruments is saved in the folder (aod_daily)[https://github.com/LudCano/aod_plotting/tree/main/aod_daily]

(*) _Important:_ this plots (or data) won't be generated (or retrieved) until the GOES satellite AWS's service is up again

Written by [Ludving Cano](mailto:lcano@chacaltaya.edu.bo) based on the code from ARSET.

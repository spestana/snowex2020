### Analysis of SnowEx 2020 snow surface temperature datasets 

The 2020 NASA SnowEx field campaign intensive observation period took place in early February 2020 at Grand Mesa in western Colorado. As part of the field campaign, a suite of snow surface temperature observations were made at a range of spatial and temporal scales to evaluate infrared imagery from GOES-R ABI. This repository contains jupyter notebooks used to analyze some of these datasets, looking particularly at observations from two cloud-free daytime periods on 8 and 11 February 2020.

### Data

All data are publicly available from the following sources:

| Dataset | Source |
| --- | --- |
| Snow Temperature Profile Time Series| [NSIDC](https://nsidc.org/data/snex20_vpts_raw/versions/1) |
| Snow Pit Measurements | [snowexsql](https://snowexsql.readthedocs.io/en/latest/) or [NSIDC](https://nsidc.org/data/snex20_gm_sp/versions/1) |
| Airborne TIR Imagery | [NSIDC](https://nsidc.org/) |
| ASTER L1T | [USGS LPDAAC](https://lpdaac.usgs.gov/products/ast_l1tv003/) via [NASA EarthData](https://www.earthdata.nasa.gov/) |
| GOES-R ABI | [NOAA GOES on AWS](https://registry.opendata.aws/noaa-goes/) via [goespy](https://github.com/palexandremello/goes-py)|

### Other resources

 * Some of this work was started as a [project](https://github.com/snowex-hackweek/hot-pow) at the 2021 [SnowEx Hackweek](https://snowex.hackweek.io/)
 * [Demo notebooks](https://github.com/spestana/snowex2020-snow-temp) for snow temperature data processing from the SnowEx 2020 field campaign

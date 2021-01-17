from pysolar.solar import get_altitude
from datetime import timezone, timedelta
import pandas as pd


def add_zenith(lat, lon, input, out=None, utc_offset=0):
    """ Add a zenith to the given dataframe"""
    df = pd.read_excel(input, na_values=-9999.0)
    df["Date/Time"] = df["Date/Time"].dt.tz_localize(timezone(timedelta(hours=utc_offset)))

    def zenith(time):
        return min(90 - get_altitude(lat, lon, time), 90)

    df["zenith"] = df["Date/Time"].apply(zenith)

    if out is None:
        out = input + "with_zenith.csv"

    df.to_csv(out, index=False)


if __name__=="__main__":
    add_zenith(lat = 51.099, lon = 10.4335784,
               input="../data/FLX_DE-Hai_FLUXNET2015_FULLSET_HH_2016-2018_beta-3_for_class.xlsx",
               out= "data/fluxnet_hainich_with_zenith.csv",
               utc_offset=1)

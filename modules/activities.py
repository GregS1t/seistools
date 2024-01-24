class Activity():
    from obspy.core import UTCDateTime
    """
    Class to create activity object to add then to ASPIC_object
    
    """
    name = ""
    starttime = UTCDateTime()
    endtime = UTCDateTime()
    sol = ""
    desc = ""
    
    def __init__(self, name="New Activity", starttime=UTCDateTime().now(), endtime=UTCDateTime().now(), desc=""):
        self.name = name
        self.starttime = starttime
        self.endtime = endtime
        self.desc = desc

    def __str__(self):
        sentence = "{} | sol: {} | starttime: {} | endtime: {} | description: {}".format(self.name, self.sol, self.starttime, self.endtime, self.desc)
        
        return sentence
        
def generate_activities_list(seis_file, sol, marslandingprop):
    """
    This function generate a list of activities for a given sol
    @seis_file: XLSX file maintained by CNES. Light green events are flagged "SEIS" and Dark green are flagged "VBB"
            This file has been manually flagged
    @sol: sol number for which we want to get activities 
    @return: A list of Activity objects
    
    """
    #
    # Open SEIS activity file to connect data with activity  
    ########################################################################

    from modules.MarsConverter import MarsConverter
    import pandas as pd
    
    def convert_string(val):
        """
        Convert Val to be useable by MarsConverter lib
        @val: SOL - float number
        """
        val = int(val)
        new_val = '{:04d}'.format(val)
        return new_val

    def manage_nan(val):
        """
        Convert Val to be useable by MarsConverter lib
        @val: 
        """
        if "nan" in val or "NaN" in val:
            val = "NaN"
        return val

    def lmst2utc(lmst_date):
        if lmst_date != "NaN":
            return marslandingprop.get_lmst_to_utc(lmst_date=lmst_date)
        else:
            return "NaN"


    seisev = None
    try: 
        fullfile = pd.read_excel(seis_file)

    except FileNotFoundError:
        print("{} not found".format(seis_file))
        
    except IOError:
        print("IO Error while reading {}".format(seis_file))
    finally:
        print("Activity file {} opened".format(seis_file))
        seisev = fullfile.loc[fullfile["Flag ASPIC"].isin(["SEIS","VBB"]), ("sol", "start LMST", "end LMST", "date UTC", "Activity", "Flag ASPIC")]
        seisev['StartLMST']= seisev['sol']
        seisev['EndLMST']= seisev['sol']

        # Concat the SOL number and the startdate -> To get a good LMST format
        seisev["StartLMST"].fillna("9999", inplace = True)
        seisev['StartLMST'] = seisev['StartLMST'].map(convert_string)+"T"+seisev['start LMST'].astype("str")+":000000"
        seisev['StartLMST'] = seisev['StartLMST'].apply(manage_nan)
        # Concat the SOL number and the enddate -> To get a good LMST format
        seisev["EndLMST"].fillna("9999", inplace = True)
        seisev['EndLMST'] = seisev['EndLMST'].map(convert_string)+"T"+seisev['end LMST'].astype("str")+":000000"
        seisev['EndLMST'] = seisev['EndLMST'].apply(manage_nan)    

        # Convert LMST Input to UTC Date
        seisev['StartUTC'] = seisev['StartLMST'].apply(lmst2utc)
        seisev['EndUTC'] = seisev['EndLMST'].apply(lmst2utc)

        print("Select activities for sol {}.".format(sol))
        seis_activity_sol = seisev.loc[(seisev['sol'] == sol)]
        seis_activity_list_sol = []
        for index in range(0, seis_activity_sol["sol"].count()):

            #print(seis_activity_sol["sol"].iloc[index], seis_activity_sol["StartUTC"].iloc[index], seis_activity_sol["EndUTC"].iloc[index], \
            #      seis_activity_sol["Activity"].iloc[index])
            activity = Activity(name="SOL_{}_{}".format(sol, index), starttime=seis_activity_sol["StartUTC"].iloc[index], endtime=seis_activity_sol["EndUTC"].iloc[index], \
                            desc=seis_activity_sol["Activity"].iloc[index])
            activity.sol = sol
            seis_activity_list_sol.append(activity)

        return seis_activity_list_sol
        
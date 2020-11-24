#!/usr/bin/env python3

# Compare Sit, Delta, and Network repo's in Python to find irregularities.

import pandas as pd  # Import the pandas package and rename it as pd.

SitInv = pd.read_csv('~/git/sit/data/inventory/inventory.csv')
DeltaRec = pd.read_csv('~/git/delta/install/receivers.csv')
Network = pd.read_csv('~/git/network/data/devices.csv')
DeltaDL = pd.read_csv('~/git/delta/install/dataloggers.csv')
DeltaRecord = pd.read_csv('~/git/delta/install/recorders.csv')

print("\nInstrument options: NetR9, Q330S/6, Obsidian\n")

Instrument = input("Enter Instrument: ")

if Instrument == "NetR9":

    # Pull NetR9's from Sit Inventory with selected columns.

	InvNetR9 = SitInv.loc[SitInv.MODEL=='NetR9',['SERIAL_NUMBER','ASSET_NUMBER','STATUS','LOCATION']]

    # Pull all NetR9's from the receiver.csv and all currently installed receivers in Delta
    # then merge using an inner join to produce a df of currently installed NetR9's.
	
	DeltaNetR9 = DeltaRec.loc[DeltaRec.Model=='TRIMBLE NETR9',['Serial','Mark']]
	DeltaRec9999 = DeltaRec.loc[DeltaRec['End Date']=='9999-01-01T00:00:00Z',['Serial','Mark']]
	Installed_NetR9 = pd.merge(DeltaNetR9, DeltaRec9999, how='inner')

    # Pull long locality names (Hostname) for NetR9's from the devices.csv in the Network repo and then merge with 
    # the currently installed NetR9's using matching site codes and and outer join to retain all columns from both df's.

	Network_NetR9 = Network.loc[Network.Model=='Trimble NetR9',['Hostname','Code']]
	NetR9_Final = pd.merge(Installed_NetR9, Network_NetR9, how='outer', left_on=['Mark'], right_on=['Code'])

    #Print Network localities that do not have a match in Delta.

	nan_rows = NetR9_Final[NetR9_Final['Serial'].isnull()]
	print("\nThe following Network repo entries are not listed in the Delta repo\n", nan_rows['Hostname'])

    #Print Delta localities that do not have a match in Network.	

	nan_rows2 = NetR9_Final[NetR9_Final['Hostname'].isnull()]
	print("\nThe following Delta repo entries are not listed in the Network repo\n", nan_rows2['Mark'])

    # Combine Sit Inventory NetR9 dataframe with previously matched Delta/Network NetR9 dataframe using common
    # serial numbers via an inner join.	

	Inv_match = pd.merge(InvNetR9, NetR9_Final, how='inner', left_on=['SERIAL_NUMBER'], right_on=['Serial'])

    # Clean up the location names from the Sit Inventory for matching the Delta/Network localities. 

	Inv_match.LOCATION = Inv_match.LOCATION.apply(
		                                          lambda x: x.lower()).apply(
		                                          lambda x: x.replace('gps','')).apply(
		                                          lambda x: x.replace('station','')).apply(
		                                          lambda x: x.replace('road','')).apply(
		                                          lambda x: x.replace('Road','')).apply(
		                                          lambda x: x.replace(" ",""))

	# Produce a true statement if Sit Inventory locality matches the Delta/Network locality and a false statement if
	# it does not. Print the false entries. 

	Inv_match['Boolean'] = Inv_match.apply(lambda x: x.LOCATION in x.Hostname, axis=1)
	A = Inv_match[~Inv_match["Boolean"]]
	print("\nThe following entries have discrepancies between the Sit Inventory and the Delta/Network Repos\n")
	with pd.option_context('display.max_rows', None, 'display.max_columns', None):
		print(A)

if Instrument == "Q330S/6":

	InvQ330S_6 = SitInv.loc[SitInv.MODEL=='Q330S/6',['SERIAL_NUMBER','ASSET_NUMBER','STATUS','LOCATION']]

	DeltaQ330S_6 = DeltaDL.loc[DeltaDL.Model=='Q330S/6',['Serial','Place']]
	DeltaDL9999 = DeltaDL.loc[DeltaDL['End Date']=='9999-01-01T00:00:00Z',['Serial','Place']]
	Installed_Q330S_6 = pd.merge(DeltaQ330S_6, DeltaDL9999, how='inner')

	Inv_match2 = pd.merge(InvQ330S_6, Installed_Q330S_6, how='inner', left_on=['SERIAL_NUMBER'], right_on=['Serial'])

	Inv_match2.LOCATION = Inv_match2.LOCATION.apply(
		                                            lambda x: x.lower()).apply(
		                                            lambda x: x.replace('gps','')).apply(
		                                            lambda x: x.replace('station','')).apply(
		                                            lambda x: x.replace('road','')).apply(
		                                            lambda x: x.replace('Road','')).apply(
		                                            lambda x: x.replace('the',''))

	Inv_match2.Place = Inv_match2.Place.apply(
		                                      lambda x: x.lower()).apply(
		                                      lambda x: x.replace('the',''))
	
	Stringedit = Inv_match2.applymap(str)

	Stringedit['cut'] = Stringedit.LOCATION.str.split().str.get(0)

	Stringedit['Boolean'] = Stringedit.apply(lambda x: x.cut in x.Place, axis=1)
	B = Stringedit[~Stringedit["Boolean"]]
	print("\nThe following entries have discrepancies between the Sit Inventory and the Delta Repo\n")
	with pd.option_context('display.max_rows', None, 'display.max_columns', None):
		print(B)

if Instrument == "Obsidian":

	InvObsidian = SitInv.loc[SitInv.MODEL=='Obsidian',['SERIAL_NUMBER','ASSET_NUMBER','STATUS','LOCATION']]

	DeltaObsidian = DeltaRecord.loc[DeltaRecord.Datalogger=='OBSIDIAN',['Serial','Station']]
	DeltaObs9999 = DeltaRecord.loc[DeltaRecord['End Date']=='9999-01-01T00:00:00Z',['Serial','Station']]
	Installed_Obsidian = pd.merge(DeltaObsidian, DeltaObs9999, how='inner')

	Network_Obs = Network.loc[Network.Model=='Kinemetrics Obsidian',['Hostname','Code']]
	Obs_Final = pd.merge(Installed_Obsidian, Network_Obs, how='outer', left_on=['Station'], right_on=['Code'])

	nan_rows3 = Obs_Final[Obs_Final['Serial'].isnull()]
	print("\nThe following Network repo entries are not listed in the Delta repo\n", nan_rows3['Hostname'])

	nan_rows4 = Obs_Final[Obs_Final['Hostname'].isnull()]
	print("\nThe following Delta repo entries are not listed in the Network repo\n", nan_rows4['Station'])

	Inv_match3 = pd.merge(InvObsidian, Obs_Final, how='inner', left_on=['SERIAL_NUMBER'], right_on=['Serial'])

	Inv_match3.LOCATION = Inv_match3.LOCATION.apply(
		                                            lambda x: x.lower()).apply(
		                                            lambda x: x.replace('gps','')).apply(
		                                            lambda x: x.replace('station','')).apply(
		                                            lambda x: x.replace('road','')).apply(
		                                            lambda x: x.replace('Road','')).apply(
		                                            lambda x: x.replace('the','')).apply(
		                                            lambda x: x.replace(" ",""))


	Inv_match3['Boolean'] = Inv_match3.apply(lambda x: x.LOCATION in x.Hostname, axis=1)
	C = Inv_match3[~Inv_match3["Boolean"]]
	print("\nThe following entries have discrepancies between the Sit Inventory and the Delta/Network Repos\n")
	with pd.option_context('display.max_rows', None, 'display.max_columns', None):
		print(C)

     









	

	



		



	
	







	

	
	

			
	
			

		
		









		




	

	



	

	






			









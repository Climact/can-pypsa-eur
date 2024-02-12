The following Excels summarize information from the optimization runs, at different granularity levels :
	- At country level / at system level
	- As temporal profiles / as capacities installed / as costs per technology

Units used 
	- installed capacities are expressed in GW_e if the technology consumes or produces electricity. An exception is made for Haber-Bosch (GW_LHV,nh3). Gas and hydrogen pipelines are on the other hand expressed in GW_LHV,ch4 and GW_LHV,h2 .
	- generation profiles are expressed in GW_e if the technology consumes or produces electricity. 
	- load profiles are expressed in MW
	- annual sums are expressed in TWh
	- costs are expressed in Euro

__________________________________________________________

The files that are available are the followings : 

- cost_countries.csv
	Cost per technology and per type of costs for all countries (EU27 and non-member states).
		* Capital costs are annualized capital costs of the asset (annualized investment costs + FOM, depending on the asset installed capacity) ;
		* Marginal costs are O&M costs of the asset (depending on the energy supplied by the asset) ;
		* Fuel costs are the extraction/import costs of the resource in the country.
	NB : transmission costs are not considered by country to avoid double counting 
___

- gas_network_capacity.csv
	Gas pipelines transmission capacities at the system level (EU27 and non-member states).
___

- gas_network_capacity_countries.csv
	Gas pipelines pipelines interconnection capacities by country for all countries (EU27 and non-member states).
___

- gas_phase_out.csv
	Installed capacities of gas powerplants > 1GW per country at EU27 level.
___

- generation_profiles.csv
	Generation profile per country, technology and timestep (timestep = 3h) over the full time horizon, expressed in MW.
___

- grid_capacity.csv
	AC and DC transmission capacities at the system level (EU27 and non-member states).
___

- grid_capacitiy_countries.csv
	AC and DC interconnection capacities by country for all countries (EU27 and non-member states).
___

- H2_network_capacity.csv
	H2 pipelines and H2-retrofitted pipelines transmission capacities at the system level (EU27 and non-member  states).
___

- H2_network_capacity_countries.csv
	H2 pipelines and H2-retrofitted pipelines interconnection capacities by country for all countries (EU27 and non-member states).
___

- import_exports.csv
	Imports and exports per countries with other countries. The values are given by year and per carrier (elec, gas and hydrogen). For each row, the exports are given per country, per year and per carrier ; imports can be read for each country by reading the values by column

Eg : NL exports in 2030 10 TWh of electricity to Belgium, 20TWh to Germany, etc

	year 	| carrier | BE 	| NL 	| DE 	| FR 	| ...
NL 	2030 	  elec	   10	  0	  20 	   0 
___

- load_profiles.csv
	Load profile per country, energy carrier (electricity, H2, NH3) and timestep (timestep = 3h) over the full time horizon, expressed in MW. An annual sum is also provided in TWh.
___

- res_potentials.csv
	RES potential per country, technology and year, expressed in GW.
___

- unit_capacities_countries.csv
	Capacities per country for all countries (EU27 and non-member states)
	NB : 'Gas', 'biomass' and 'uranium' technologies represent the extraction/import capacities of gas, biomass and uranium. 
___

- h2_production.png
	Profile of h2 production, storage and use (for fuel-cells) for all countries per timestep (timestep = 3h) over the full time horizon and over the month of January.

___

- storage_units.png
	Main short-term and long-term storage units state of charge for all countries per timestep (timestep = 3h) over the full time horizon for long-term storage and over the month of January for short term storage
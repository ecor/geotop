# GEOtop river network codes, pixel borders definition keywords and discharge output variables description

## Pixel border set up in GEOtop
```
FreeDrainageAtBottom=0
```		  
FreeDrainageAtBottom=0 means that at the bottom there is an impermable condifion (no flux drains down and no water has been lost from the GEOtop 3D computational domain), otherwise FreeDrainageatBottom=1 means water drains down for gravity and is lost.
```
FreeDrainageAtLateralBorder=1
```		
FreeDrainageAtLateralBorder=1 means  water can flow by gravity laterally following the topography  gradient with no gradient of soil water pressure head, there is no impermeable barrier at the lateral border of the catchment soil.  By the pratical point of view, in case of mountain divide, this condition could not affect the results significantly, beacause soil water piezometric head should have a maximum (gradient equal to zero) , respect to the impermeable lateral boundary condition (FreeDrainageAtBottom=0). 

Which cells are included in the calculation, as our DEM has a buffer of NaNs to not have a shift with the slope and aspect file…are  the NaNs included?

Some tips:

* If you start with an already cutted DEM and no other maps of derived topographic properties are present (i.e. slope, drainage direction ...), then GEOtop calculates them my himself. However, in this case it turns to novalue every cell close to the border of your DEM , not only the first / last rows of the dem matrix. This could be taken in account, especially for small catchments.

## Channels network definitions in GEOtop

There are several codes to define land or channel network pixels.


* **Code -1 => LAND pixel**: an incoming discharge from outside is considered (as rain/irrigation). Be aware that this option has been little used and tested.

* **Code 0 => LAND pixel**: this is the typical definitio for an hillslope pixel. If this pixel is on the border:
	- the border is impermeable
	- no lateral surface runoff 

* **Code 1 => LAND pixel**. If this pixel is on the border:
	- only lateral surface runoff is allowed

* **Code 2 => LAND pixel** If this pixel is on the border it behaves weir-wise ("stramazzo"):
	- the border is permeable above an user-defined elevation in the saturated part, defined as the keyword "DepthFreeSurfaceAtTheBoundary" [mm]. Be aware that this setting has been little tested.

* **10 => CHANNEL pixel** This is the usual code for for the river network. If it is on the border it becomes the river network outlet:
	- water can exit only in the river network (Qtot)
	- the border is impermeable outside the river.

* **11 => CHANNEL pixel** If it is on the border it behaves both as channel pixel and as pixel_type 1 (river discharge and lateral surface runoff)

* **12 => CHANNEL pixel** If it is on the border it behaves oth as channel pixel and as pixel_type 2 (river discharge and lateral surface and subsurface runoff)

You can see how we specified the river network
close to basin outlet in the borden example:
https://github.com/geotopmodel/geotop/blob/master/tests/3D/Borden05m/input-maps/net05.asc


## Meaning of the different discharge terms in the file "DischargeFile"
 
1. **DATE[day/month/yearhour:min]**
2. **t[days]** Days from simulation begin
3. **JDfrom0** Time as MATLAB, serial date number (starting from 1 is January 1, 0000)
4. **JD** Julian Day
5. **Qtot[m3/s]** Total discharge from the river network outlet(s)
6. **Vsup/Dt[m3/s]** total surface water flux exchanging the river network with the surrounding soil in the considered time step
7. **Vsub/Dt[m3/s]** total sub-surface water flux exchanging the river network with the surrounding soil
8. **Vchannel[m3]** total water volume present in the river network
9. **Qoutlandsup[m3/s]** surface runoff at the catchment lateral borders outside the river network.
10. **Qoutlandsub[m3/s]** sub-surface runoff at the catchment lateral borders outside the river network.
11. **Qoutbottom[m3/s]** total water flux percolating at the bottom boundary of the domain (i.e. deep groundwater recharge).

 


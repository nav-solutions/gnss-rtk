GNSS-RTK
========

[![crates.io](https://img.shields.io/crates/v/gnss-rtk.svg)](https://crates.io/crates/gnss-rtk)
[![Rust](https://github.com/rtk-rs/gnss-rtk/actions/workflows/rust.yml/badge.svg)](https://github.com/rtk-rs/gnss-rtk/actions/workflows/rust.yml)
[![crates.io](https://docs.rs/gnss-rtk/badge.svg)](https://docs.rs/gnss-rtk)
[![rustc](https://img.shields.io/badge/rustc-1.64%2B-blue.svg)](https://img.shields.io/badge/rustc-1.64%2B-blue.svg)

[![License](https://img.shields.io/badge/license-MPL_2.0-orange?style=for-the-badge&logo=mozilla)](https://github.com/rtk-rs/qc-traits/blob/main/LICENSE)

Position Velocity Time (PVT) solution solver in Rust.

GNSS-RTK is an easy and efficient navigation solver that supports both PPP and RTK navigation. 
This makes it the ideal solution for most navigation applications.  
It is easy to deploy, the configuration is simple and even the default setup will exhibit good results
(no tweaking involved!).

GNSS-RTK is flexible and efficient:

* you can navigate with a single signal in sight
* you don't have to sample L1 and can navigate with modern signals
* it supports navigation without phase range
* a special CPP method for dual frequency pseudo range (no phase range)
which behaves like a slow converging PPP method
* is a true surveying tool because it can operate without apriori knowledge
* it can fulfill the challenging task of RTK / Geodetic reference station calibration
by deploying a complete PPP survey

Other cool features:
  * works in all supported timescales
  * can navigate using a conic azimuth mask (min and max azimuth angle).
  In this case, we only select vehicles from that very region of the compass.
  * could potentially apply to other Planets, if we make some function more generic
  and propose better atmosphere interfaces.

GNSS-RTK does not care about SV or signal modulations. It cares
about physics, distances, frequencies and environmental phenomena.
This means you operate it from whatever data source you have at your disposal,
as long as you can provide the required inputs.

Applications
============

GNSS-RTK is used by the following applications

* [Post processed PPP with RINEX-Cli (app)](https://github.com/georust/rinex) 
* [Post processed Common View Timing with RINEX-Cli (app)](https://github.com/georust/rinex) 
* [Real Time PPP with RT-NAVI (app)](https://github.com/rtk-rs/rt-navi)

Examples and more information
=============================

* Some examples are shipped with this repo, they will teach you how to deploy the solver,
at least in basic setups
* [The RINEX Wiki](https://github.com/georust/rinex/wiki) describes extensive application of this framework, at a high level

Framework
=========

GNSS-RTK includes itself within and is closely tied to the following libraries:

* [ANISE](https://github.com/nyx-space/anise) for orbital calculations
* [Nyx-space](https://github.com/nyx-space/nyx) for advanced navigation
* [Hifitime](https://github.com/nyx-space/hifitime) for timing
* [GNSS-rs](https://github.com/rtk-rs/gnss) for basic GNSS definitions
* Nalgebra for all calculations

Licensing
=========

This library is part of the [RTK-rs framework](https://github.com/rtk-rs) which
is delivered under the [Mozilla V2 Public](https://www.mozilla.org/en-US/MPL/2.0) license.

Work in Progress
================

This framework is in active development. Most features have been stabilized
and exhibit very good performances, some are still unavailable or needs to be improved.

The following topics needs to be addressed (also, refer to the active Github Issues):

* `ppp` phase navigation is not fully completed yet. 
We achieve very good performance by performing long survey. 
PPP will allow much faster convergence.
* `rtk-ppp` is not fully possible yet, only basic RTK exists. 

Refer to this portal (Wiki pages, Discussions..) and the RINEX Wiki to understand what this tool is capable of.

PVT Solutions
=============

The Solver's objective is to resolve precise PVT solutions.
A minimum of 4 SV must be observed in standard navigation.  
When navigating in Fixed Altitude mode, only 3 SV must be observed.   
When navigating in Time Only mode, a single SV needs to be observed.   

When performing a survey (read dedicated paragraph), 4 SV needs to be observed until
the solver fully initializes itself. Use the returned object (PVTSolution or Error) to determine
whether you can relax the constraint on sky observations (yet) or not.

The preset criteria are manually set in the configuration file (or config script).
At the moment, refer to the RINEX Wiki or RINEX scripts database, for meaningful examples.

Depending on the preset configuration, other requirements will apply to the previous list, most importantly:

- `CPP` strategy will required pseudo range observation on a secondary frequency
- `PPP` strategy will required pseudo range and phase observations on two frequencies
- SNR, Elevation and Azimuth mask will require to gather the required amount of SV within those conditions

Each PVT solution contains the Dilution of Precision (DOP) and other meaningful information, like which SV
contributed to the solution. We have the capability to express the clock offset in all supported Timescale.

Strategy and other settings
===========================

The solver's behavior and output results are highly dependent on the [selected strategy](https://docs.rs/gnss-rtk/latest/gnss_rtk/prelude/enum.Method.html).

Advanced strategies require deeper knowledge and most likely more tuning of the solver configuration. 
The Rust/JSON infrastructure is powerful enough though, to allow to only define the [config parts](https://docs.rs/gnss-rtk/latest/gnss_rtk/prelude/struct.Config.html)
you are interested in: others will simply default.

[PVTSolutionType](https://docs.rs/gnss-rtk/latest/gnss_rtk/prelude/enum.PVTSolutionType.html) defines the type of solutions we want to form and therefore,
the minimum amount of SV we need to gather. As previously stated, other criteria like `min_sv_elev` or `max_sv_azim` will restrict the condition on those vehicles that they must fit in
to be considered. 

When `fixed_altitude` is set to a certain value, the quantity of required SV is reduced by 1.  
This has no impact when `PVTSolutionType` is set to `TimeOnly`.

The `SolverOpts` configuration gives more advanced options on how to tweak the solver. Briefly, this allows to

- select one of our Navigation Filters, like Kalman filter or LSQ
- define the PVT solutions confirmation criteria

`Modeling` defines what physical and environmental phenomena we compensate for.   
Modeling are closely tied to the selected solver strategy. For example, 
models that impact at the centimetric level like the sunlight rate, are not meaningful in strategies other than advanced PPP.
On the other hand, you will not reach metric solutions, whatever the strategy might be, if a minimum of physical phenomena are not accounted for.

Atmosphere and Environmental biases
===================================

This solver is always capable of modeling all conditions and form a solution.   
It is important to understand how our API is designed and operate it the best you can to get the best results.

Troposphere Delay
==================

Troposphere bias always needs to be estimated.  
By default, the solver will use a model implemented in the [model::tropo API].  
If you're in position to determine yourself the Tropospherical Delay components (TropoComponents structure)
at the required latitude and _Epoch_, you are highly encouraged to provide your data.
To do so, we use a function pointer that can serve as a TropoComponents source.
TropoComponents evaluation parameters (function pointer arguments) should be :

- _Epoch_
- altitude (above sea level) expressed in meters
- latitude expressed as decimal degrees

For _Epochs_ where the data source is not capable to supply data, that is not a problem, we will rely on the internal model.

Example of handmade TropoComponents provider :  

TODO

:warning: When using the internal model, we recommend applying a >= 5° elevation mask to the candidates.

Ionosphere Delay
================

TODO

A priori position
=================

The solver can be armed with a priori knowledge (rough idea of the final position),
or can operate in complete autonomy. In this case, the solver will initialize itself
very accurately, this requires one extra step.

Getting started
===============

Refer to our example applications to understand how to operate our API in more detail.

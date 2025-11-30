GNSS-RTK
========

[![Rust](https://github.com/nav-solutions/gnss-rtk/actions/workflows/rust.yml/badge.svg)](https://github.com/nav-solutions/gnss-rtk/actions/workflows/rust.yml)
[![Rust](https://github.com/nav-solutions/gnss-rtk/actions/workflows/daily.yml/badge.svg)](https://github.com/nav-solutions/gnss-rtk/actions/workflows/daily.yml)  

[![crates.io](https://img.shields.io/crates/v/gnss-rtk.svg)](https://crates.io/crates/gnss-rtk)
[![crates.io](https://docs.rs/gnss-rtk/badge.svg)](https://docs.rs/gnss-rtk)
[![crates.io](https://img.shields.io/crates/d/gnss-rtk.svg)](https://crates.io/crates/gnss-rtk)
[![discord server](https://img.shields.io/discord/1342922474110586910?logo=discord)](https://discord.gg/EqhEBXBmJh)

[![MRSV](https://img.shields.io/badge/MSRV-1.82.0-orange?style=for-the-badge)](https://github.com/rust-lang/rust/releases/tag/1.82.0)
[![License](https://img.shields.io/badge/license-AGPLv3.0-orange?style=for-the-badge)](https://github.com/nav-solutions/gnss-rtk/blob/main/LICENSE)

The `GNSS-RTK` library provides Position Velocity Time (PVT) solution solvers,
with abstract and flexible interfaces that may apply to most navigation scenarios
and applications. Whether they are real-time or post-processed does not matter.

This library is thoroughly validated, including the accuracy of the solutions
are validated for each release.

<div align="center">
    <p>
        <br>Static CPP (Galileo, E1+E5b)</br> surveying a professional geodetic marker
    </p>
    <a href=https://github.com/nav-solutions/rinex-cli/blob/main/plots/front-page/map.png>
        <img src=https://github.com/nav-solutions/rinex-cli/blob/main/plots/front-page/map.png alt="Plot">
    </a>
</div>

<div align="center">
    <p>
        Static CPP (Galileo, E1+E5b) residual errors between GNSS-RTK
and profesionnaly calibrated position.
    </p>
    <a href=https://github.com/nav-solutions/rinex-cli/blob/main/plots/front-page/coordinates.png>
        <img src=https://github.com/nav-solutions/rinex-cli/blob/main/plots/front-page/coordinates.png alt="Plot">
    </a>
</div>

<div align="center">
    <p>
        <br>CPP (Galileo E1+E5b)</br>
pedestrian profile (no reference on the ground)
    </p>
    <a href=https://github.com/nav-solutions/rinex-cli/blob/main/plots/front-page/roaming-ppp1.png>
        <img src=https://github.com/nav-solutions/rinex-cli/blob/main/plots/front-page/roaming-ppp1.png alt="Plot">
    </a>
</div>

Licensing
=========

This library is released under [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.en.html).

GNSS-RTK / P.V.T solutions
==========================

The main objective of this library is to propose efficient GNSS solvers, that
perform the rather complex task of solving Position, Velocity and Timing (P.V.T) solutions,
possibly from scratch (without initial guess).

Indeed, GNSS navigation allows to obtain the position of the receiver in 3D (in space) 
and possibly in time as well.

The timing solutions (4th component) is the state of the receiver in the navigation timescale,
and is obtained by PPP navigation technique (not RTK). Our timing solutions
may be expressed in all supported Timescales and we offer an interface for precise timing corrections.

For the spatial coordinates, whether the receiver is moving or not is application dependent and one
of the major point of the navigation filter is to track the target in all cases. To help the navigation filter,
it is recommended to tune your `UserParameters` according to your application. This includes

- `UserProfile` for the type of application (basically, the expected average velocity).
Selecting an incorrect profile will not prohibit the filter to converge to the truth, it just reduces
the convergence time.

- `ClockProfile`

The proposed API is compatible with both real-time and post-processed navigation, so it should be usable
in any navigation application. For real-time applications, we offer interesting features:

- the user profile is updated at eatch `Epoch`, if you can monitor your clock profile and have
erratic application, you can take advantage of that.

- when the RTK network goes down, you can switch to PPP technique temporarily

GNSS-RTK is currently limited to ground navigation on planet Earth. Although it should be feasible to
make it more abstract and suitable for navigation on other planets, it is not scheduled as of today.
If you want to this happen, feel free to open discussions.

Summary
=======

Select a navigation method, depending on available signals and desired accuracy:

| Method        | Physics                                  | Accuracy      |  Application                                            |
|---------------|------------------------------------------|---------------|---------------------------------------------------------|
| `SPP`         | Single Frequency Pseudo Range navigation |  2/5          | Low cost devices, Degraded setups                       |
| `CPP`         | Dual Frequency Pseudo Range navigation   |  4/5          | Mid cost dual frequency systems, Timing applications    |
| `PPP`         | Dual Frequency Pseudo Range + Phase      |  5/5          | High cost devices, Timing applications, Profesionnal surveying, RTK Stations calibration.. | 

For each set of synchronous measurements (oftentimes referred to as "Observations"), you can then either 

1. apply the RTK Navigation technique (Differential navigation, relying on a single ground reference),
by using `Solver::rtk()`.

2. or PPP Navigation technique (Absolute navigation without reference on the ground),
by using `Solver::ppp()`.

RTK should be prefered for geometric applications, especially roaming applications.
So far, we experience between x2-x5 improvements of the geometric accuracy when comparing absolute
to differential runs. In RTK, the most important thing is to maintain a short baseline (distance between
the rover and the reference).

Our RTK API will not resolve the clock state in any case. If you are interested in the timing
solutions, you should prefer PPP runs. When the RTK network goes down and the reference
goes unreachable, you can safely switch to `Solver::ppp` temporarily, the navigation filter continues
and iterates from the past state at all times.

Summary of supported modes:

| Principle  | API               | Method | Latest version | Accuracy |
|------------|-------------------|--------|----------------|----------|
| `RTK`      | `Solver::rtk_run` | SPP    |  OK            |  2/5     |
|------------|-------------------|--------|----------------|----------|
| `RTK`      | `Solver::rtk_run` | CPP    |  OK            |  3/5     |
|------------|-------------------|--------|----------------|----------|
| `RTK`      | `Solver::rtk_run` | PPP    |  OK            |  5/5     |
|------------|-------------------|--------|----------------|----------|
| `PPP`      | `Solver::ppp_run` | SPP    |  OK            |  1/5     |
|------------|-------------------|--------|----------------|----------|
| `PPP`      | `Solver::ppp_run` | CPP    |  OK            |  2/5     |
|------------|-------------------|--------|----------------|----------|
| `PPP`      | `Solver::ppp_run` | PPP    |  NOK           |  4/5     |
|------------|-------------------|--------|----------------|----------|

NOK: unreleased, unavailable or untested. Will at best diverge, and most likely panic if you deploy this mode.

Solver, strategies and API
==========================

The current version proposes a unique object called `Solver` that is compatible with both absolute
and differential (1 reference) navigation. Once your measurements have been collected,
you select the navigation technique you want to deploy with either `ppp_run` and `rtk_run`.

Whether the rover is moving or not does not matter, you should simply tune your `UserParameters`,
accordingly. But because RTK improves the geometric accuracy, it will give best results for moving rovers.

Our Method structure is there to define the signals you want to use. This allows us to
give total freedom to the end user while allow both absolute and differential navigation.
Absolute navigation is not particularly tied to phase signal, those are two different things.
And you can use our solver to do RTK without phase observations as well.

Because this API is real-time oriented, you provide measurements for each Epoch.
Assuming people are mostly interested by highest accuracy, you would naturally select 
RTK+PPP or PPP. But you can catch errors and adapt. Let's say L5 becomes unavailable for some time,
you can continue navigating but temporarily degrade the selected Method to SPP.

Application Programming Interface (API)
=======================================

The API requires you implement of several `Traits` to provide information
that the solving process demands:

- `OrbitSource` must be implemented to provide the orbital states of each satellites
at a specifing point in time. It is the only mandatory traits, all the others can be tied to zero
and you will still obtain a (poor) solution

- `EnvironmentalBias` must be implemented to describe the environmental perturbations. 
On earth ground navigation, it comprises the ionospheric and tropospheric perturbation.
Depending on the selected mode, some of it may be disregarded. Starting from `CPP` mode,
the ionospheric delay is physically canceled: if you tie this bias to zero, you will still obtain
a solution of good quality. For the tropospheric delay, you can either implement one of the models we propose
or implement your own. It is mandatory to do so to obtain a high quality solution;

- `SpacebornBias` must be implemented to obtain a basic solution. Indeed, if you don't provide the
`SatelliteClockCorrection` at the requested point in time, the solution will be off by thousands of kilometers.

Each traits are wrapped inside an `Rc<>` (Reference Counter) which makes this API compatible with real-time navigation.
This allows your structure(s) to live elsewhere, with so called "internal multability" (read on this pseudo concept)
that the data collection process demands.

RTK-Trait
=========

Each Ground station must implement the `RTKBase` trait, to describe its current test.  
For each set of rover measurements, we request synchronous remote observations. 
In practice, especially in real-time applications, it is impossible for the observations
to truly be synchronous. The important aspect is to reduce the time difference (also referred to as
the RTK "aging") to a minimum. Considering the ground reference is expected to remain static,
this will work just fine. For moving base, this becomes vital. In summary, the latest state
you described for the ground reference should be as close to the truth as you can, and any error
will reflect on the accuracy of the rover solution. 

We request for a "name" for each base station, this makes the logs more meaningful
and is introduced in advanced even though we do not allow more than 1 reference per Epoch.

GNSS-RTK is currently limited to a single ground reference.

This library has not been tested with a moving base yet, but if you can describe your
base correctly at all times, you should be able to obtain correct results.

Orbit Provider
==============

When measurements have been proposed, the solver will issue a request to this function pointer,
soon after. That you must answer, otherwise this candidate proposal will be dropped and will not contribute
to the solution. That is not a big deal, assuming you are in position to provide many measurements.
But we will always need at least 4 proposal to pass this step and future "post-fit" attitude filters
that you may have selected.

Other function pointers are much easier to implement and are not mandatory.

Environment models & perturbations
==================================

We provide a function pointer that allows you to apply the Ionosphere and Troposphere perturbations
compensation. If you are in position to apply your own database and/or model, you then have an interface
to do so. Returning no compensation will simply degrade the quality of your solution.

If you don't have a custom model, simply pick one of our internal models and apply it as the answer
to the request. For example, `TroposphereModel` can serve as a simple answer to the troposphere perturbation request.

Time transposition function
===========================

We require a function pointer to perform any time transposition. This is very useful if you have
access to correction tables for finer temporal corrections. In case you don't, simply return
the transposition using `Epoch.to_time_scale()`.

This API allows us to support all timescales supported by Hifitime and obtain higher
precision for people who have access to such data.

Examples
========

This library is no longer shipped with examples, for the simple reason that they are complicated to maintain,
and our framework now offers the two kind of applications you can think of to this library:

- [RT-Navi](https://github.com/nav-solutions/rt-navi) for real-time navigation, which is obviously hardware dependent,
and involves multi-threading.
- [GNSS-Qc](https://github.com/nav-solutions/gnss-qc) for post-processed navigation, which is more complex

Logs
====

This library uses `$RUST_LOG` to report the current state, what it is doing, analyze your preset
and give recommendations, and finally, help debug potential use cases.

For correct debugging, you should set `$RUST_LOG=debug`. For ulimtate debugging, `$RUST_LOG=trace`
will give all the information we output.

Deployment
==========

:warning: ANISE requires internet access either

1. at built time, when built with `embed_ephem` option 
2. at first deployment ever, when using low precision models and built without `embed-ephem` option
3. regularly, when using highest precision models, regardless of the compilation options.

This library exposes the `embed_ephem` compilation option. By using it, it is possible to run up to
mid-range and precise application without internet access. Ultra precise applications require
regular updates from the Internet.

Configuration simplicity
========================

Despite solving complex stuff, we strive to offer simplistic interfaces and configuration setups,
which is a challenge. Our `Config` structure is fully `serde` compatible,
and is easy to deploy and customize, refer to the online documentation:

Notes on Navigation Method
==========================

Absolute or differential is selected at deployment. The 
[navigation method](https://docs.rs/gnss-rtk/latest/gnss_rtk/prelude/enum.Method.html) is defined
by the configuration preset.

- `SPP`: is a Pseudo Range navigation method for single signal / degraded / limited setups.
Phase range observations are disregarded, unless you intend to use the L1+C1 enhancing smoothing technique.
You can enhance the `SPP` with an external model for environment perturbations. You can use the models we provide 
to answer that requirement. If you do not provide an accurate model (which is a complex thing to do), you cannot obtain
accurate solutions. We recommend using a different technique for any setups that allow you to do so.
It is often times interesting to degrade a CPP or PPP run back to SPP to actually see the huge difference
these two made.

- `CPP`: starting at CPP, we use physical cancellation of the ionosphere bias.
You need to provide two separate frequencies for this mode to operate.
CPP is like SPP and purely based on Pseudo Range navigation, phase range are disregarded,
unless you intend to use the L1+C1+L2+C2 smoothing technique.

- `PPP`: has the same requirements and objectives, but phase observations
are now also expected. It can be further enhanced by deploying the code smoothing technique.
Since PPP is very heavy and requires all signals to be available, our PPP presets activate the 
smoothing technique by default.

All of these are independent of the rover's profile AND the navigation technique:

- you can use any signal strategy along any rover profile
- you an use any signal strategy with either static or dynamic applications
- you can use any signal strategy along RTK differential technique

Enhanced Navigation techniques
==============================

We provide a few options to enhance the accuracy of the solutions, depending of the context and input scenario.

- When using `SPP` technique, you can provide L1 phase observations and deploy the L1+C1 code smoothing technique,
to improve the accuracy of your results. This is very interesting for degraded / low-cost use cases.

- When using the `CPP` technique, you can provide L1+L2 phase observations and deploy the L1+C1+L2+C2 code smoothing
technique, which does not really make sense, because you are providing all the necessary input to switch to PPP technique.
You should switch to `PPP` technique any time L1 + L2 phase observations are feasible

- Enhancing the `PPP` technique with code smoothing will improve the accuracy of the solutions.

Surveying and Apriori Knowledge
===============================

`GNSS-RTK` is a complete surveying tool that can operate without apriori knowledge.  
In otherwords, it is possible to obtain a precise solution from scratch.  
[Refer to this page that demonstrates this capacity](https://github.com/nav-solutions/rinex-cli/blob/main/demos/SURVEY.md). 

Signal and Measurement Flexibility
==================================

One objective is to make this solver flexible and capable to adapt to most exploitation scenarios.  
Answering 100% of usecases is impossible, for the simple reason that GNSS applications
cover a very wide spectrum and vary a lot. Yet we provide central key elements that are very interesting:

* Possibility to navigate using a single signal. Obviously, the samples you provide
should match your navigation technique at all times. Yet, we have several navigation techniques,
some are compatible with single frequency observation. This makes our solver compatible with
degraded or low-cost environments

* Military codes: our library only cares about frequencies.
Whether you arrive from a decoded military or civilian signal does not matter.
L1 is treated as main reference as always.

* High precision frequencies: E6 (Galileo) and LEX (QZSS) are both known

* :warning: This library requires L1 frequency to be present for advanced CPP or PPP techniques
(as reference signal). We allow L2 or L5 as subsidary signal, without prioritizing them.

Constellations, Timecales & Absolute time
=========================================

Like other topics, `GNSS-RTK` tries to be convenient to operate
and flexible, regarding the temporal solution.

We support all timescales provided by `Nyx-Space/Hifitime`.   
The `Time` trait should be implemented for applications that require very precise temporal solutions.  
It allows following the timescale states precisely (using external corrections).  

This library is flexible enough to let you express your measurements in the timescale you want,
and express the solution in the timescale you want. UTC timescale is fully supported.

Other features
==============

Non exhaustive list of other interesting features

- It is possible to apply a custom conic Azimutal + Elevation mask

Framework
=========

`GNSS-RTK` relies on the following libraries:

* Nalgebra
* [ANISE](https://github.com/nyx-space/anise) 
* [Hifitime](https://github.com/nyx-space/hifitime) 

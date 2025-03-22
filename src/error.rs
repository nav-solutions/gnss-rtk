use thiserror::Error;

use anise::{
    almanac::{metaload::MetaAlmanacError, planetary::PlanetaryDataError},
    errors::{AlmanacError, PhysicsError},
};

use crate::prelude::{Epoch, SV};

#[derive(Debug, PartialEq, Error)]
pub enum Error {
    /// Not enough candidates were proposed, with respect to navigation parameters.
    #[error("not enough candidates provided")]
    NotEnoughCandidates,
    /// Survey initialization (no apriori = internal guess)
    /// requires at least 4 SV in sight temporarily, whatever
    /// your navigation technique.
    #[error("survey initialization requires at least 4 SV temporarily")]
    NotEnoughInitializationCandidates,
    /// PreFit (signal quality, other..) criterias
    /// have been applied but we're left with not enough vehicles that match
    /// the navigation technique: no attempt.
    #[error("not enough candidates match pre-fit criteria")]
    NotEnoughPreFitCandidates,
    /// PostFit (state solver and other) have been resolved,
    /// but we're left with not enough vehicles that match
    /// the navigation technique: no attempt.
    #[error("not enough candidates match post-fit criteria")]
    NotEnoughPostFitCandidates,
    /// Failed to parse navigation method
    #[error("non supported/invalid strategy")]
    InvalidStrategy,
    #[error("not enough post-fit candidates to form a matrix")]
    MatrixMinimalDimension,
    #[error("internal error: invalid matrix setup")]
    MatrixDimension,
    #[error("failed to form matrix (invalid input or not enough data)")]
    MatrixFormationError,
    /// Invalid orbital states or bad signal data may cause the algebric calculations
    /// to wind up here.
    #[error("failed to invert matrix")]
    MatrixInversion,
    /// Invalid orbital states or bad signal data may cause the algebric calculations
    /// to wind up here.
    #[error("resolved time is `nan` (invalid value(s))")]
    TimeIsNan,
    /// Invalid orbital states or bad signal data may cause the algebric calculations
    /// to abort.
    #[error("internal navigation error")]
    NavigationError,
    /// Failed to initialize navigation filter
    #[error("nav filter initialization error")]
    NavigationFilterInitError,
    #[error("missing pseudo range observation")]
    MissingPseudoRange,
    /// [Method::CPP] requires the special signal combination to exist.
    /// This require the user to sample PR on two separate frequencies.
    #[error("failed to form pseudo range combination")]
    PseudoRangeCombination,
    /// [Method::PPP] requires the special signal combination to exist.
    /// This require the user to sample PR + PH on two separate frequencies.
    #[error("failed to form phase range combination")]
    PhaseRangeCombination,
    /// Each [Candidate] state needs to be resolved to contribute to any PPP resolution attempt.
    #[error("unresolved candidate state")]
    UnresolvedState,
    /// Each [Candidate] presented to the Bancroft solver needs a resolved state.
    #[error("bancroft requires 4 fully resolved candidates")]
    UnresolvedStateBancroft,
    /// When [Modeling.sv_clock_bias] is turned on and we attempt PPP resolution,
    /// it is mandatory for the user to provide [ClockCorrection].
    #[error("missing clock correction")]
    UnknownClockCorrection,
    /// Physical non sense due to bad signal data or invalid orbital state, will cause us
    /// abort with this message.
    #[error("physical non sense: rx prior tx")]
    PhysicalNonSenseRxPriorTx,
    /// Physical non sense due to bad signal data or invalid orbital state, will cause us
    /// abort with this message.
    #[error("physical non sense: t_rx is too late")]
    PhysicalNonSenseRxTooLate,
    // /// Solutions may be invalidated and are rejected with [InvalidationCause].
    // #[error("invalidated solution, cause: {0}")]
    // InvalidatedSolution(InvalidationCause),
    /// In pure PPP survey (no RTK, no position apriori knowledge = worst case scenario),
    /// [Solver] is initiliazed by [Bancroft] algorithm, which requires
    /// temporary 4x4 navigation and pseudo range sampling (whatever your navigation technique),
    /// until at least initialization is achieved.
    #[error("bancroft solver error: invalid input ?")]
    BancroftError,
    /// [Bancroft] initialization process (see [BancroftError]) will wind up here
    /// in case unrealistic or bad signal observation or orbital states were forwarded.
    #[error("bancroft solver error: invalid input (imaginary solution)")]
    BancroftImaginarySolution,
    /// PPP navigation technique requires phase ambiguity to be solved prior any attempt.
    /// It is Okay to wind up here for a few iterations, until the ambiguities are fixed
    /// and we may proceed to precise navigation. We will reject solving attempt until then.
    /// Hardware and external events may reset the ambiguity fixes and it is okay to need to
    /// rerun through this phase for a short period of time. Normally not too often, when good
    /// equipment is properly operated.
    #[error("unresolved signal ambiguity")]
    UnresolvedAmbiguity,
    /// [Solver] requires [Almanac] determination at build up and may wind-up here this step is in failure.
    #[error("issue with Almanac: {0}")]
    Almanac(AlmanacError),
    /// [Solver] uses local [Almanac] storage for efficient deployments
    #[error("almanac setup issue: {0}")]
    MetaAlmanac(MetaAlmanacError),
    /// [Solver] requires to determine a [Frame] from [Almanac] and we wind-up here if this step is in failure.
    #[error("frame model error: {0}")]
    EarthFrame(PlanetaryDataError),
    /// Any physical non sense detected by ANISE will cause us to abort with this error.
    #[error("physics issue: {0}")]
    Physics(PhysicsError),
    /// Remote observation is required for a [Candidate] to contribute in RTK solving attempt.
    /// You need up to four of them to resolve. We may print this internal message and still
    /// proceed to resolve, as [SV] may go out of sight of rover or reference site.
    #[error("missing observation on remote site {0}({1})")]
    MissingRemoteRTKObservation(Epoch, SV),
    /// In RTK resolution attempt, you need to observe all pending [SV] on reference site as well.
    /// If that is not the case, we abort with this error.
    #[error("missing observations on remote site")]
    MissingRemoteRTKObservations,
    #[error("unknown or non supported frequency")]
    InvalidFrequency,
    #[error("rejected troposhere delay: model divergence?")]
    RejectedTropoDelay,
    #[error("rejected ionosphere delay: model divergence?")]
    RejectedIonoDelay,
    #[error("nav filter converged to physicaly invalid state")]
    StateUpdate,
}

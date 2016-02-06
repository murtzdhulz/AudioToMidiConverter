% Intelligent Sound Processing Toolbox

% Demonstrations
%  ISP_AUDIOREADDEMO      - Demonstration of how to use isp_audioread
%  ISP_DISTANCEMEASUREDEMO - Demonstrate various aspects of music distance measures
%  isp_distributedemo     - Demonstration of isp_distribute.
%  ISP_EVALUATEDEMO1      - Demonstration of how to use the evaluation framework
%  ISP_EVALUATEDEMO2      - Demonstration of how to use the evaluation framework
%  ISP_GMMDEMO            - Demonstration that uses the Gaussian mixture model functions.
%  ISP_MIDIDEMO           - Demonstration of midi input/output
%  ISP_RHYTHMDEMO         - Demonstration of isp_tirhythm and the evaluation framework

% Classification
%  ISP_GMMDEMO            - Demonstration that uses the Gaussian mixture model functions.
%  ISP_GMMDISTANCE        - Compute distance between two Gaussian mixture models.
%  ISP_GMMPROBABILITY     - Probability of vectors according to a specified Gaussian Mixture Model.
%  ISP_GMMRAND            - Generate random vectors from a Gaussian Mixture Model.
%  ISP_GMMTRAIN           - Train a Gaussian mixture model from data
%  ISP_KMEANS             - Binary tree k-means algorithm

% Distributed Computing
%  ISP_DISTRIBUTEDEMO_HELPER1 - Helper function for isp_distributedemo
%  ISP_DISTRIBUTEDEMO_HELPER2 - Helper function for isp_distributedemo
%  isp_distributedemo     - Demonstration of isp_distribute.
%  ISP_DISTRIBUTE         - Distribute job for parallel execution on multiple machines.
%  ISP_RUNJOB             - Execute pending, distributed jobs.

% Evaluation
%  ISP_COMPUTEDISTANCE    - Compute distance matrix.
%  isp_evaluatedistancematrix - Helper function for isp_evaluate
%  ISP_EVALUATE           - Evaluate a musical distance measures
%  isp_extractfeature     - Extract features from a song.
%  ISP_MAKESONGLIST       - Generate a song list from various collections
%  ISP_PLOTRESULTS        - Plots results returned by isp_evaluate
%  ISP_TESTSKELETON       - Internal function providing a framework for isp_evaluate.

% Features
%  ISP_AR                 - Compute multivariate autoregressive coefficients.
%  ISP_ARWRAPPER          - Wrapper function for isp_ar.
%  ISP_CHROMA             - Compute chroma vectors
%  ISP_IFCHROMAGRAM       - Compute chromagram from instantaneous frequency
%  ISP_IFGRAM             - Instantaneous frequency by phase derivative.
%  ISP_IFPTRACK           - Pitch track based on instantaneous frequency.
%  ISP_MELFILTERBANK      - Compute the mel scaled filterbank.
%  ISP_MFCCAT             - Reimplementation of MFCC from the Auditory Toolbox
%  ISP_MFCCGMMKL          - Definition of MFCC-GMM-Kullback-Leibler distance measure
%  ISP_MFCCSIG            - Compute the MFCC of a signal.
%  ISP_MFCCVB             - MFCC implementation based on Mike Brookes' Voicebox
%  ISP_SPECTROGRAM        - Compute the spectrogram of a signal.
%  ISP_TICHROMA           - Define time scale invariant chroma-based distance measure.
%  ISP_TIRHYTHM           - Define time scale insensitive measure of rhythmic distance.

% Sampled audio
%  ISP_AUDIOREADDEMO      - Demonstration of how to use isp_audioread
%  ISP_AUDIOREAD          - Read audio file.
%  ISP_MP3READ            - Read wave data and song information from an MP3 file.
%  ISP_WAVMODIFY          - Modifies properties of a WAV struct.

% MIDI
%  ISP_MIDIDEMO           - Demonstration of midi input/output
%  ISP_MIDIGMNAME         - Return names of instruments in General MIDI Level 1 sound set.
%  ISP_MIDIMODIFY         - Modifies properties of a MIDI struct.
%  ISP_MIDIREAD           - Read MIDI file
%  ISP_MIDISHOW           - Visualize a midi song
%  ISP_MIDISYNTH          - Synthesize a MIDI song to wave data.
%  ISP_MIDIWRITE          - Write MIDI file.

% Miscellaneous
%  ISP_ASSERT             - Print error message if assertion does not hold.
%  ISP_BSEARCH            - Binary search
%  ISP_CALLEXECUTABLE     - Locate and execute external command (internal ISP toolbox command).
%  ISP_INTERPRETARGUMENTS - Interpret arguments and set default values for unspecified fields.
%  ISP_LINESTYLE          - Create options for line and marker styles.
%  ISP_NONEMPTY           - Return the first non-empty input argument
%  ISP_PICKLE             - Convert Matlab data to column vector
%  ISP_PLOT               - Plot many-dimensional variable
%  ISP_SHOWSPECTROGRAM    - Plots a spectrogram
%  ISP_SUBPLOTLABEL       - Plot per-figure labels
%  ISP_TEMPFILE           - Generate temporary file name
%  ISP_TOOLBOXPATH        - Returns the path to the ISP Toolbox
%  ISP_UNPICKLE           - Recover original Matlab object from pickled column vector.

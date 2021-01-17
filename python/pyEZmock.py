#!/usr/bin/env python3
import os
from warnings import warn

class pyEZmock:

  def __init__(self, workdir=None,
      exe='/global/u2/z/zhaoc/work/pyEZmock/bin/EZmock.sh',
      pk_exe='/global/u2/z/zhaoc/work/pyEZmock/bin/POWSPEC.sh',
      xi_exe='/global/u2/z/zhaoc/work/pyEZmock/bin/FCFC_2PT_BOX.sh',
      bk_exe='/global/u2/z/zhaoc/work/pyEZmock/bin/BISPEC_BOX.sh'):
    if workdir is None:
      raise ValueError('Working directory should be set via `workdir`.')
    if not (os.path.isfile(exe) and os.access(exe, os.X_OK)):
      raise IOError(f'Invalid EZmock executable: {exe}')
    if not (os.path.isfile(pk_exe) and os.access(pk_exe, os.X_OK)):
      raise IOError(f'Invalid power spectrum executable: {exe}')
    if not (os.path.isfile(xi_exe) and os.access(xi_exe, os.X_OK)):
      raise IOError(f'Invalid 2PCF executable: {exe}')
    if not (os.path.isfile(bk_exe) and os.access(bk_exe, os.X_OK)):
      raise IOError(f'Invalid bispectrum executable: {exe}')

    if not os.path.isdir(workdir):
      try: os.mkdir(workdir)
      except: raise IOError(f'Cannot create working directory: {workdir}')

    self.workdir = workdir
    self.exe = exe
    self.pk_exe = pk_exe
    self.xi_exe = xi_exe
    self.bk_exe = bk_exe

    # Initialise EZmock parameters
    self.param = dict(
      boxsize = None,
      num_grid = None,
      redshift = None,
      num_tracer = None,
      pdf_base = None,
      dens_scat = None,
      rand_motion = None,
      dens_cut = None,
      seed = 1,
      omega_m = 0.307115,
      init_pk = '/global/u2/z/zhaoc/work/pyEZmock/data/PlanckDM.linear.pk'
    )

    # Initialise clustering settings
    self.pk = False             # indicate whether to compute power spectra
    self.pkl = []               # power spectrum multipoles to be evaluated
    self.pkgrid = 256           # grid size for power spectra evaluation
    self.kmax = 0.3             # maximum k for the power spectra
    self.dk = 0.01              # bin size of k for the power spectra
    self.pk_ref = None          # reference power spectrum
    self.pk_col = []            # columns of the reference P(k) multipoles

    self.xi = False             # indicate whether to compute 2PCFs
    self.xi_z = []              # 2PCF multipoles to be evaluated
    self.rmax = 150             # maximum r for the 2PCFs
    self.dr = 5                 # bin size of r for the 2PCFs
    self.xi_ref = None          # reference 2PCF
    self.xi_col = []            # columns of the reference 2PCF multipoles

    self.bk = False             # indicate whether to compute bispectrum
    self.bkgrid = 256           # grid size for bispectra evaluation
    self.k1 = [None, None]      # range of the k1 vector for the bispectrum
    self.k2 = [None, None]      # range of the k2 vector for the bispectrum
    self.bkbin = 20             # number of output bins for the bispectrum
    self.bk_ref = None          # reference bispectrum
    self.bk_col = None          # column of the reference bispectrum

    # Other parameters
    self.history_dict = []      # history of evaluated parameter sets
    self.odir = None            # output directory for the current run
    self.script = None          # script for the current run
    self.ezfile = None          # EZmock catalogue for the current run


  def set_param(self, boxsize, num_grid, redshift, num_tracer,
      pdf_base=None, dens_scat=None, rand_motion=None, dens_cut=None,
      seed=1, omega_m=0.307115,
      init_pk='/global/u2/z/zhaoc/work/pyEZmock/data/PlanckDM.linear.pk'):
    """
    Set parameters for EZmock evaluation.

    Parameters
    ----------
    boxsize: float
        Side length of the cubic periodic box.
    num_grid: int
        Number of grids per side for the density field.
    redshift: float
        Final redshift of the catalogue.
    num_tracer: int
        Number of tracers to be generated.
    pdf_base: float, optional
        Base number for PDF mapping.
    dens_scat: float, optional
        Density modification parameter.
    rand_motion: float, optional
        Parameter for the random local motion.
    dens_cut: float, optional
        The critical density.
    seed: int, optional
        Random seed.
    omega_m: float, optional
        Density parameter at z = 0.
    init_pk: str, optional
        Initial power spectrum.
    Reference
    ---------
    https://arxiv.org/abs/2007.08997
    """
    self.param['boxsize'] = float(boxsize)
    self.param['num_grid'] = int(num_grid)
    self.param['redshift'] = float(redshift)
    self.param['num_tracer'] = int(num_tracer)
    if not pdf_base is None: self.param['pdf_base'] = float(pdf_base)
    if not dens_scat is None: self.param['dens_scat'] = float(dens_scat)
    if not rand_motion is None: self.param['rand_motion'] = float(rand_motion)
    if not dens_cut is None: self.param['dens_cut'] = float(dens_cut)
    self.param['seed'] = int(seed)
    if (not init_pk is None) and (not os.path.isfile(init_pk)):
      raise IOError(f'init_pk does not exist: {init_pk}')
    self.param['init_pk'] = init_pk
    self.param['omega_m'] = float(omega_m)
    if self.param['omega_m'] <= 0 or self.param['omega_m'] > 1:
      raise ValueError('omega_m must be between 0 and 1')


  def set_clustering(self, pk=False, pk_ell=[0,2], pk_grid=256, pk_kmax=0.3,
      pk_dk=0.01, pk_ref=None, pk_ref_col=[5,6], xi=False, xi_ell=[0,2],
      xi_rmax=150, xi_dr=5, xi_ref=None, xi_ref_col=[4,5], bk=False,
      bk_grid=256, bk_k1=[0.04,0.06], bk_k2=[0.09,0.11], bk_nbin=20,
      bk_ref=None, bk_ref_col=4):
    """
    Set configurations for clustering measurements.

    Parameters
    ----------
    pk: bool, optional
        Indicate whether to compute power spectrum.
    pk_ell: tuple of ints, optinal
        Legendre multipoles of power spectrum to be evaluated.
    pk_grid: int, optional
        Grid size for power spectra evaluations.
    pk_kmax: float, optional
        Maximum k for the power spectra.
    pk_dk: float, optional
        Bin size of k for the power spectra.
    pk_ref: str, optional
        File for the reference power spectrum. The first column must be k.
    pk_ref_col: tuple of ints, optional
        Columns (counting from 0) of the reference power spectrum multipoles.

    xi: bool, optional
        Indicate whether to compute 2-point correlation function (2PCF).
    xi_ell: tuple of ints, optinal
        Legendre multipoles of 2PCF to be evaluated.
    xi_rmax: float, optional
        Maximum separation for the 2PCFs.
    xi_dr: float, optional
        Bin size of separation for the 2PCFs.
    xi_ref: str, optional
        File for the reference 2PCF. The first column must be r.
    xi_ref_col: tuple of ints, optional
        Columns (counting from 0) of the reference 2PCF multipoles.

    bk: bool, optional
        Indicate whether to compute bispectrum.
    bk_grid: int, optional
        Grid size for bispectra evaluations.
    bk_k1: (float, float), optional
        Range of the first k vector for the bispectra.
    bk_k2: (float, float), optional
        Range of the second k vector for the bispectra.
    bk_nbin: int, optional
        Number of output bins for the bispectra.
    bk_ref: str, optional
        File for the reference bispectrum. The first column must be theta.
    bk_ref_col: int, optional
        Column of the reference bispectrum.
    """
    self.pk = pk
    if self.pk and len(pk_ell) == 0: raise ValueError('pk_ell is empty')
    self.pkl = pk_ell
    self.pkgrid = pk_grid
    self.kmax = pk_kmax
    self.dk = pk_dk
    if (not pk_ref is None) and (not os.path.isfile(pk_ref)):
      raise IOError(f'pk_ref does not exist: {pk_ref}')
    self.pk_ref = pk_ref
    if not pk_ref is None and len(pk_ell) != len(pk_ref_col):
      raise ValueError('pk_ell and pk_ref_col should have equal length')
    self.pk_col = pk_ref_col

    self.xi = xi
    if self.xi and len(xi_ell) == 0: raise ValueError('xi_ell is empty')
    self.xil = xi_ell
    self.rmax = xi_rmax
    self.dr = xi_dr
    if (not xi_ref is None) and (not os.path.isfile(xi_ref)):
      raise IOError(f'xi_ref does not exist: {xi_ref}')
    self.xi_ref = xi_ref
    if not xi_ref is None and len(xi_ell) != len(xi_ref_col):
      raise ValueError('xi_ell and xi_ref_col should have equal length')
    self.xi_col = xi_ref_col

    self.bk = bk
    self.bkgrid = bk_grid
    self.k1 = bk_k1
    self.k2 = bk_k2
    self.bkbin = bk_nbin
    if (not bk_ref is None) and (not os.path.isfile(bk_ref)):
      raise IOError(f'bk_ref does not ebkst: {bk_ref}')
    self.bk_ref = bk_ref
    self.bk_col = int(bk_ref_col)


  def run(self, nthreads, queue=None, walltime=30, partition='haswell',
      boxsize=None, num_grid=None, redshift=None, num_tracer=None,
      pdf_base=None, dens_scat=None, rand_motion=None, dens_cut=None,
      seed=None, omega_m=None, init_pk=None):
    """
    Run the job for EZmock generation and clustering measurements.

    Parameters
    ----------
    nthreads: int
        Number of OpenMP threads used for the run.
    queue: str, optional
        Queue of the job to be submitted (e.g. 'debug' and 'regular').
        If not provided, the job script has to be run manually.
    walltime: int, optional
        Limit on the total run time (in minutes) of the jobs.
    partition: str, optional
        Specify the architecture of the nodes for the job.
        It has to be 'haswell' or 'knl'.
    The rest of the parameters are the same as those for `set_param`.
    """
    import copy
    from subprocess import Popen, PIPE

    previous_param = copy.deepcopy(self.param)

    if not boxsize is None: self.param['boxsize'] = float(boxsize)
    if not num_grid is None: self.param['num_grid'] = int(num_grid)
    if not redshift is None: self.param['redshift'] = float(redshift)
    if not num_tracer is None: self.param['num_tracer'] = int(num_tracer)
    if not pdf_base is None: self.param['pdf_base'] = float(pdf_base)
    if not dens_scat is None: self.param['dens_scat'] = float(dens_scat)
    if not rand_motion is None: self.param['rand_motion'] = float(rand_motion)
    if not dens_cut is None: self.param['dens_cut'] = float(dens_cut)
    if not seed is None: self.param['seed'] = int(seed)
    if not omega_m is None:
      self.param['omega_m'] = float(omega_m)
      if self.param['omega_m'] <= 0 or self.param['omega_m'] > 1:
        raise ValueError('omega_m must be between 0 and 1')
    if not init_pk is None:
      if not os.path.isfile(init_pk):
        raise IOError(f'init_pk does not exist: {init_pk}')
      self.param['init_pk'] = init_pk

    if None in self.param.values():
      raise ValueError('Please set EZmock parameters via `set_param` or `run`.')

    # Create the path and name of the job script for the current run
    bname = self.__get_bname(self.param)
    self.odir = '{}/{}'.format(self.workdir, bname)
    if not os.path.isdir(self.odir):
      try: os.mkdir(self.odir)
      except: raise IOError(f'Cannot create directory: {self.odir}')
    self.script = '{}/run_{}.sh'.format(self.odir, bname)
    if os.path.isfile(self.script):
      raise IOError('The job script exists, please wait if the job ' \
          'has already been submitted, or submit/run the script manually: \n'
          f'{self.script}')
    self.ezfile = f'EZmock_{bname}.dat'

    # Store the previous set of parameters as history
    if not None in previous_param.values():
      prev_bname = self.__get_bname(previous_param)
      if os.path.isfile(f'{self.workdir}/{prev_bname}/DONE'):
        self.history_dict.append(previous_param)

    # Generate contents of the job script file.
    nthreads = int(nthreads)
    if nthreads <= 0: raise ValueError(f'Invalid nthreads: {nthreads:d}')
    jobstr = ('#!/bin/bash\n#SBATCH -n 1\n#SBATCH -L SCRATCH\n'
        f'#SBATCH -o {self.odir}/stdout_%j.txt\n'
        f'#SBATCH -e {self.odir}/stderr_%j.txt\n')
    if not queue is None:
      jobstr += f'#SBATCH -q {queue}\n'
      if partition != 'haswell' and partition != 'knl':
        raise ValueError('Invalid job partition: {partition}')
      jobstr += f'#SBATCH -C {partition}\n#SBATCH -c {nthreads}\n'
    if not walltime is None:
      jobstr += '#SBATCH -t {:d}\n'.format(int(walltime))

    jobstr += f'\nexport OMP_NUM_THREADS={nthreads}\n\ncd {self.odir}\n\n'

    run_mock = True
    if os.path.isfile(f'{self.odir}/{self.ezfile}'):
      warn(f'EZmock will not be run, as file exists: {self.odir}/{self.ezfile}')
      run_mock = False
    if run_mock: jobstr += self.__mock_cmd(bname)

    if self.pk: jobstr += self.__pk_cmd()
    if self.xi: jobstr += self.__xi_cmd()
    if self.bk: jobstr += self.__bk_cmd()

    jobstr += f'echo 1 > DONE\n'

    # Save the job script and submit it if applicable
    with open(self.script, 'w') as f: f.write(jobstr)

    if queue is None:
      print('Job script generated. Please run the following command manually:')
      print(f' \nbash {self.script}')
    else:
      process = Popen(['/usr/bin/sbatch',self.script], shell=False,
          stdout=PIPE, stderr=PIPE, text=True)
      sts = process.wait()
      for line in process.stdout: print(line, end='')
      for line in process.stderr: print(line, end='')
      if sts != 0:
        print('Job submission failed. Please resubmit the script manually:')
        print(f'sbatch {self.script}')
        print(f'or run the script directly:\nbash {self.script}')


  def plot(self):
    """
    Plot the clustering measurements of the previous runs, the references,
    and the current run.
    """
    import numpy as np
    import matplotlib.pyplot as plt

    # Determine the number of subplots
    nplot = 0
    if self.pk: nplot += len(self.pkl)
    if self.xi: nplot += len(self.xil)
    if self.bk: nplot += 1
    if nplot == 0:
      print('No clustering measurements specified.')
      print('Please set via `set_clustering`')
      return

    # Create subplots
    ncol = 3
    if nplot <= 4: ncol = nplot
    nrow = int(np.ceil(nplot / ncol))
    figw = min(15, nplot * 5)
    figh = 3 * nrow
    fig = plt.figure(figsize=(figw,figh))
    plt.subplots_adjust(wspace=0.3, hspace=0.3)
    ax = [plt.subplot2grid((nrow,ncol), (i//ncol,i%ncol)) for i in range(nplot)]
    for a in ax: a.grid(ls=':', c='dimgray', alpha=0.6)

    iplot = 0
    alpha_hist = 0.8
    ls_ref = 'k:'
    ls_curr = 'k-'

    # Plot power spectra multiples
    if self.pk:
      # Plot histories first
      for j, hist in enumerate(self.history_dict):
        bname = self.__get_bname(hist)
        ifile = f'{self.workdir}/{bname}/PK_EZmock_{bname}.dat'
        d = np.loadtxt(ifile, unpack=True)
        for i, ell in enumerate(self.pkl):
          ax[iplot+i].plot(d[0], d[5+i]*d[0]**1.5, alpha=alpha_hist,
              label=f'history {j:d}')
      # Plot the reference
      if not self.pk_ref is None:
        d = np.loadtxt(self.pk_ref, unpack=True)
        for i, c in enumerate(self.pk_col):
          ax[iplot+i].plot(d[0], d[c]*d[0]**1.5, ls_ref, label='ref')
      # Plot the current results
      ifile = f'{self.odir}/PK_{self.ezfile}'
      if os.path.isfile(ifile):
        d = np.loadtxt(ifile, unpack=True)
        for i, ell in enumerate(self.pkl):
          ax[iplot+i].plot(d[0], d[5+i]*d[0]**1.5, ls_curr, label='current')
      else: warn('the current job may not have been finished')
      # Set axis labels and ranges
      for i, ell in enumerate(self.pkl):
        ax[iplot+i].set_xlabel(r'$k$')
        ax[iplot+i].set_ylabel(r'$k^{{1.5}} P_{:d} (k)$'.format(ell))
        ax[iplot+i].set_xlim(0, self.kmax)
      iplot += len(self.pkl)

    # Plot 2PCF
    if self.xi:
      # Plot histories first
      for j, hist in enumerate(self.history_dict):
        bname = self.__get_bname(hist)
        ifile = f'{self.workdir}/{bname}/2PCF_EZmock_{bname}.dat'
        d = np.loadtxt(ifile, unpack=True)
        for i, ell in enumerate(self.xil):
          ax[iplot+i].plot(d[0], d[3+i]*d[0]**2, alpha=alpha_hist,
              label=f'history {j:d}')
      # Plot the reference
      if not self.xi_ref is None:
        d = np.loadtxt(self.xi_ref, unpack=True)
        for i, c in enumerate(self.xi_col):
          ax[iplot+i].plot(d[0], d[c]*d[0]**2, ls_ref, label='ref')
      # Plot the current results
      ifile = f'{self.odir}/2PCF_{self.ezfile}'
      if os.path.isfile(ifile):
        d = np.loadtxt(ifile, unpack=True)
        for i, ell in enumerate(self.xil):
          ax[iplot+i].plot(d[0], d[3+i]*d[0]**2, ls_curr, label='current')
      else: warn('the current job may not have been finished')
      # Set axis labels and ranges
      for i, ell in enumerate(self.xil):
        ax[iplot+i].set_xlabel(r'$s$')
        ax[iplot+i].set_ylabel(r'$s^{{2}} \xi_{:d} (s)$'.format(ell))
        ax[iplot+i].set_xlim(0, self.rmax)
      iplot += len(self.xil)

    if self.bk:
      # Plot histories first
      for j, hist in enumerate(self.history_dict):
        bname = self.__get_bname(hist)
        ifile = f'{self.workdir}/{bname}/BK_EZmock_{bname}.dat'
        d = np.loadtxt(ifile, unpack=True)
        ax[iplot].plot(d[0]/np.pi, d[4], alpha=alpha_hist,
            label=f'history {j:d}')
      # Plot the reference
      if not self.bk_ref is None:
        d = np.loadtxt(self.bk_ref, unpack=True)
        ax[iplot].plot(d[0]/np.pi, d[4], ls_ref, label='ref')
      # Plot the current results
      ifile = f'{self.odir}/BK_{self.ezfile}'
      if os.path.isfile(ifile):
        d = np.loadtxt(ifile, unpack=True)
        ax[iplot].plot(d[0]/np.pi, d[4], ls_curr, label='current')
      else: warn('the current job may not have been finished')
      # Set axis labels
      ax[iplot].set_xlabel(r'$\theta_{12} / \pi$')
      ax[iplot].set_ylabel(r'$B (\theta_{12})$')
      iplot += 1

    ax[-1].legend()


  def params(self):
    """
    Print the current set of EZmock parameters.
    """
    for key, value in self.param.items():
      print(key, '=', value)


  def history(self):
    """
    Print the histories of the EZmock parameter sets.
    """
    for i, param in enumerate(self.history_dict):
      print(f'{i:d}:', param)


  def clear(self, slicer):
    """
    Clear history entries defined by `slicer`

    Parameters
    ----------
    slicer:
        The slice of the histories to be cleared.
        It must be generated using the `slice` function.
    """
    if type(slicer).__name__ != 'slice':
      raise TypeError('slicer must be generated using the `slice` function')
    del(self.history_dict[slicer])


  def __get_bname(self, param):
    """
    Generate the basename of files for a given parameter set.

    Parameters
    ----------
    param: dict
        Dictionary storing a set of EZmock parameters.

    Return
    ------
    The basename as a string.
    """
    if None in param.values(): return None
    bname = (f"B{param['boxsize']:g}"
        f"G{param['num_grid']:d}"
        f"Z{param['redshift']:g}"
        f"N{param['num_tracer']:d}_"
        f"b{param['pdf_base']:g}"
        f"d{param['dens_scat']:g}"
        f"r{param['rand_motion']:g}"
        f"c{param['dens_cut']:g}")
    return bname


  def __mock_cmd(self, bname):
    """
    Generate the command for generating the redshift space EZmock.

    Parameters
    ----------
    bname: str
        Basename of the EZmock catalogue.

    Returns
    -------
    The command as a string
    """
    from cosmoEZ import flatLCDM

    # Compute structure growth parameters
    cosmo = flatLCDM(omega_m = self.param['omega_m'])
    z1 = 1 + self.param['redshift']
    a = 1. / z1
    grow2z0 = cosmo.growth2z0(a)
    hubble = cosmo.hubble(a)
    zdist = cosmo.growthf(a) * hubble * a

    # Generate the command for running EZmock
    jobstr = ("echo '&EZmock_v0_input\n"
        f"datafile_prefix = \"EZmock_{bname}_real\"\n"
        f"datafile_path = \"./\"\n"
        f"iseed = {self.param['seed']:d}\n"
        f"boxsize = {self.param['boxsize']:g}\n"
        f"grid_num = {self.param['num_grid']:d}\n"
        f"redshift = {self.param['redshift']:g}\n"
        f"grow2z0 = {grow2z0:g}\n"
        f"expect_sum_pdf = {self.param['num_tracer']:d}\n"
        f"expect_A_pdf = {self.param['pdf_base']:g}\n"
        f"density_cut = {self.param['dens_cut']:g}\n"
        f"scatter2 = {self.param['dens_scat']:g}\n"
        f"zdist_rate = {zdist:g}\n"
        f"zdist_fog = {self.param['rand_motion']:g}\n"
        "density_sat = 100\nscatter = 10\nmodify_pk = 0.0\n"
        "modify_pdf = 0\nantidamping = 2\n"
        "use_whitenoise_file = .false.\nwhitenoise_file = \"\"\n"
        f"pkfile = \"{self.param['init_pk']}\"\n"
        f"pknwfile = \"{self.param['init_pk']}\"\n"
        "compute_CF = .false.\ncompute_CF_zdist = .false.\n"
        "dilute_factor = 0.3\nskiplines = 0\ntwod_corr_suffix = \"\"\n"
        f"max_r = 50\nbin_size = 5\nom = {self.param['omega_m']:g}\n/'"
        f" | {self.exe} || exit\n\n")

    # Generate the command for applying redshift space distortions
    bsize = self.param['boxsize']
    # Check boundaries and remove non-numerical entries (such as nan)
    jobstr += ("awk '{CONVFMT=\"%.8g\";OFMT=\"%.8g\"; "
        f"if ($1>=0 && $1<{bsize:g} && $2>=0 && $2<{bsize:g}"
        f" && $3>=0 && $3<{bsize:g} && $6+0==$6) print $1,$2,"
        f"($3+$6*{z1:g}/{hubble:g}+{bsize:g})%{bsize:g} }}' "
        f"EZmock_{bname}_real.dat > {self.ezfile} || exit\n\n")

    return jobstr


  def __pk_cmd(self):
    """
    Generate the command for computing power spectrum of EZmock.

    Returns
    -------
    The command as a string
    """
    bsize = self.param['boxsize']
    poles = '[' + ','.join(f'{i:d}' for i in self.pkl) + ']'
    ifile = self.ezfile
    ofile = f'PK_{self.ezfile}'

    jobstr = (f'{self.pk_exe} -d {ifile} --data-formatter "%lf %lf %lf" '
        f"-p '[$1,$2,$3]' -s T -B {bsize:g} -G {self.pkgrid:d} -n 1 -i F "
        f"-l '{poles}' -k 0 -K {self.kmax:g} -b {self.dk:g} -a {ofile} || "
        "exit\n\n")

    return jobstr


  def __xi_cmd(self):
    """
    Generate the command for computing 2-point correlation function of EZmock.

    Returns
    -------
    The command as a string
    """
    bsize = self.param['boxsize']
    poles = '[' + ','.join(f'{i:d}' for i in self.xil) + ']'
    ifile = self.ezfile
    ofile = f'2PCF_{self.ezfile}'

    jobstr = (f"{self.xi_exe} -i {ifile} -l D -f '%lf %lf %lf' -x '[$1,$2,$3]' "
        f'-b {bsize:g} -B 1 -p DD -P {ofile}.dd -e "DD/@@-1" -E {ofile}.xi2d '
        f"-m '{poles}' -M {ofile} --s-min 0 --s-max {self.rmax:g} "
        f'--s-step {self.dr:g} --mu-num 60 --dist-prec 0 -S 0 || exit\n\n')

    return jobstr


  def __bk_cmd(self):
    """
    Generate the command for computing power spectrum of EZmock.

    Returns
    -------
    The command as a string
    """
    bsize = self.param['boxsize']
    ifile = self.ezfile
    ofile = f'BK_{self.ezfile}'

    jobstr = (f'{self.bk_exe} -i {ifile} -b 0 -B {bsize:g} -g {self.bkgrid:d} '
        f'-w 1 -x 0 -p {self.k1[0]:g} -P {self.k1[1]:g} -q {self.k2[0]:g} '
        f'-Q {self.k2[1]:g} -n {self.bkbin:d} -o {ofile} -y || exit\n\n')

    return jobstr


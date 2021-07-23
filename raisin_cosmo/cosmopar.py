import getdist.plots as gplot
import pylab as plt
import numpy as np
import pylab as plt

def getw_cosmosis(name=''):

    # should run
    # postprocess {inifile} -o plots -p {name} --no-plots
    # if needed
    with open(f"plots/{name}_means.txt") as fin:
        for line in fin:
            if line.startswith("cosmological_parameters--omega_m"):
                omegam,omegamerr = float(line.split()[1]),float(line.split()[2])
            elif line.startswith("cosmological_parameters--w"):
                w,werr = float(line.split()[1]),float(line.split()[2])
    return omegam,omegamerr,w,werr

def getw(name=''):

    g = gplot.getSinglePlotter(chain_dir='/n/holystore01/LABS/berger_lab/Lab/djones01/RAISIN/chains/chains/')
    samples = g.sampleAnalyser.samplesForRoot(name,settings={'ignore_rows':0.9})

    p = samples.getParams()

    plt.close()
    print(f"w: {samples.mean(p.w)} +/- {samples.std(p.w)}")
    print(f"Om: {samples.mean(p.omegam)} +/- {samples.std(p.omegam)}")
    return(samples.mean(p.w),samples.std(p.w),
           samples.mean(p.omegam),samples.std(p.omegam))

def geth0(name=''):

    g = gplot.getSinglePlotter(chain_dir='/scratch/midway2/rkessler/djones/cosmomc/chains_2015/chains/')
    samples = g.sampleAnalyser.samplesForRoot(name)

    p = samples.getParams()

    plt.close()
#    import pdb; pdb.set_trace()
    print(f"H0: {samples.mean(p.__dict__['H0'])} +/- {samples.std(p.__dict__['H0'])}")
    print(f"w: {samples.mean(p.w)} +/- {samples.std(p.w)}")
    print(f"Om: {samples.mean(p.omegam)} +/- {samples.std(p.omegam)}")
    return(samples.mean(p.w),samples.std(p.w),
           samples.mean(p.omegam),samples.std(p.omegam))


def getom(name=''):

    g = gplot.getSinglePlotter(chain_dir='/scratch/midway2/rkessler/djones/cosmomc/chains_2015/chains/')
    samples = g.sampleAnalyser.samplesForRoot(name)

    p = samples.getParams()

    plt.close()
    print(f"Om: {samples.mean(p.omegam)} +/- {samples.std(p.omegam)}")
    return(samples.mean(p.omegam),samples.std(p.omegam))

def cosmosys(postprocess=False):

    syslist = ['stat','all','photcal','massdivide','biascor','pecvel','kcor','mwebv','lowzcal','hstcal','tmpl','lcfitter']
    titles = ['Stat.','All Sys.','Phot. Cal.','Mass Step',
              'Bias Corr.','Pec. Vel.','$k$-corr.',
              'MW E(B-V)','Low-$z$ Cal','$HST$ Cal','Template Flux','NIR SN Model']
    fileprefix='raisin'


    if postprocess:
        os.system("postprocess RAISIN_stat.ini -o plots -p raisin_stat --no-plots")
        os.system("postprocess RAISIN_all.ini -o plots -p raisin_all --no-plots")

        for s in syslist:
            os.system(f"postprocess RAISIN_{s}.ini -o plots -p raisin_{s} --no-plots")
        
    wstat,werrstat,omstat,omerrstat = getw_cosmosis(name='raisin_stat')
    wall,werrall,omall,omerrall = getw_cosmosis(name='raisin_all')

    tblhdr = """\\begin{deluxetable}{lccc}
\\tabletypesize{\\scriptsize}
\\tablewidth{0pt}
\\tablecaption{Summary of Systematic Uncertainties on $w$}
\\tablehead{Error&$\\Delta w$&$\\sigma_w$&$\\sigma_w/\\sigma_{\\textrm{stat}}$}"""
    print(tblhdr)

    for s,t,i in zip(syslist,titles,range(len(titles))):
        f = '%s_%s'%(fileprefix,s)
        w,werr,om,omerr = getw_cosmosis(name=f)
        if werr < werrstat: 
            werrsys = 0.0
        else:
            werrsys = np.sqrt(werr**2. - werrstat**2.)
        if f == 'foundps1_all':
            print('%s&%.3f&%.3f&%.3f\\\\'%(
                t,wall-wstat,werrsys,werrsys/werrstat))
        else:
            print('%s&%.3f&%.3f&%.3f\\\\'%(
                t,w-wstat,werrsys,werrsys/werrstat))

        if i == 0:
            print('\\tableline\\\\')

        tblfooter = """\\startdata
\\enddata
\\label{table:syserr}
\\end{deluxetable}"""
    print(tblfooter)

def syspiechart(ax=None,sysval=[0.028,0.045,0.008,0.026,0.044,0.010,0.010,0.009],
                title=None,
                syslist=['Phot. Cal.','Bias Corr.','$k$-corr.',
                         'Pec. Vel.','Mass\nStep','NIR SN\nModel','MW E(B-V)','Template Flux'],
                explode=[0,0,0,0,0,0,0,0],radius=1.4,fontsize=13,makebold=False,startangle=80):
    import matplotlib.patheffects as path_effects

    systot = np.sum(sysval)
    sysval /= np.sum(sysval)
    colors = ['#fbb4ae',
              '#b3cde3',
              '#ccebc5',
              '#decbe4',
              '#fed9a6',
              '#ffffcc',
              '#e5d8bd',
              '#fddaec']

    def absolute_value(val):
        #a    = np.round(val/np.sum(sysval), 0)
        if val*systot/100 > 0.004:
            return '%.3f\n(%i%%)'%(val*systot/100,val)
        else:
            return ''
    plt.rcParams['font.size'] = 17
    if not ax:
        plt.clf()
        ax = plt.axes()
    patches, texts, autotexts = ax.pie(sysval, labels=syslist, colors=colors,
                                       autopct=absolute_value, shadow=False, startangle=startangle,
                                       labeldistance=1.08,explode=explode,
                                       wedgeprops = {'linewidth': 2, 'edgecolor':'k'},
                                       pctdistance=0.85,radius=radius)
    ax.set_title(title,y=1.025)
    for patch in patches:
        patch.set_path_effects([path_effects.SimpleLineShadow(),
                                path_effects.Normal()])
        
    for label,text in zip(autotexts,texts):
        if makebold:
            if 'CC' in text._text or 'SALT2' in text._text or 'Bias' in text._text:
                label._fontproperties._weight = 'bold'
                text._fontproperties._weight = 'bold'
        label._fontproperties._size = fontsize
    plt.savefig('raisin_syspie.png',dpi=200)

def getcorner(name=''):

    import corner

    g = gplot.getSinglePlotter(chain_dir='/scratch/midway2/rkessler/djones/cosmomc/chains_2015/chains/')
    samples = g.sampleAnalyser.samplesForRoot(name)

    p = samples.getParams()

    mat = np.zeros([len(p.w),3])
    mat[:,0] = p.__dict__['H0']
    mat[:,1] = p.w
    mat[:,2] = p.omegam
    corner.corner(mat,labels=['H$_0$','$w$','$\Omega_m$'],
				  quantiles=[0.16,0.5,0.84],show_titles=True,
				  title_kwargs={"fontsize":20},label_kwargs={"fontsize":20},bins=40)
    plt.savefig('cosmo_corner.png',dpi=200)
    #plt.close()


if __name__ == "__main__":
    getw('raisin_all_planck18')
    #getw('planck18')
    #getw('raisin_stat')
    #geth0('raisin_all')
    #getcorner('raisin_all')
    #getom('RAISIN_all_ocdm')
    #getom('RAISIN_all_lcdm')
    #getw('sn_cmb_omw_0')

    #cosmosys()
    #syspiechart()

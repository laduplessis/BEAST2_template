import os, sys, io, yaml, json, datetime
import numpy as np
from fnmatch import fnmatch
from optparse import OptionParser
from Bio import SeqIO

usage = "usage: %prog [option]"
parser = OptionParser(usage=usage)

parser.add_option("-i","--inputpath",
                  dest = "inputpath",
                  default = "",
                  metavar = "path",
                  help = "Path to input files [required]")

parser.add_option("-c","--config",
                  dest = "config",
                  default = "*.cfg",
                  metavar = "",
                  help = "Pattern to match for config files [default = %default]")

parser.add_option("-x","--template",
                  dest = "template",
                  default = "",
                  metavar = "path",
                  help = "Path to template XML file [required]")

parser.add_option("-o","--outputpath",
                  dest = "outputpath",
                  default = "",
                  metavar = "path",
                  help = "Path to save output file in [required]")

parser.add_option("-n","--name",
                  dest = "name",
                  default = "",
                  metavar = "path",
                  help = "Name of the runs [required]")

(options,args) = parser.parse_args()

if (options.inputpath != ""):
	config         = options.config
	inputpath      = os.path.abspath(options.inputpath)+"/"
else:
	config         = options.config[options.config.rfind("/")+1:]
	inputpath      = os.path.abspath(options.config[:options.config.rfind("/")])+"/"

################################################################################################################################  


def writeAlignment(filename, output, msa_id="msa"):
      """
      	Writes alignment in BEAST format to output buffer
      	Return dictionary with sequence ids and dates
      """

      datedict = dict()
      for seq in SeqIO.parse(filename, "fasta"):		

            # Process date
            seqdate = seq.id.split('_')[-1]
            #(seqdate, sequpper, seqlower)  = getYearDate(datestring, datefmt = "%d-%m-%Y")

            #print("%s\t%f\t%f\t%f\n" % (datestring, seqdate, sequpper, seqlower))

            #if (sequpper - seqlower != 0): 
            #    sys.stdout.write("Estimating sampling date for: %s\n" % seq.id)
            #    tipdates[seq.id] = [sequpper, seqlower]
            #

            datedict[seq.id] = seqdate

            # Write sequence
            output.write('\t\t\t<sequence id="%s:%s" taxon="%s" totalcount="4" value="%s"/>\n' % (seq.id, msa_id, seq.id, seq.seq))

      #

      return datedict

      #return({'datetrait' : datedict, 'tipdates' : tipdates})
#

def writeDateTrait(dates, output):
      """
      	Writes dateTrait in BEAST format for sequence ids in date dictionary
      """

      maxdate = max(dates.values())
      mindate = min(dates.values())

      traits = []
      for seqid in dates:
            traits.append('\n\t\t\t\t\t%s=%.13f' %(seqid, dates[seqid]))
            if (dates[seqid] == mindate): 
                  sys.stdout.write("Most recent sample: %s\n" % seqid)
            if (dates[seqid] == maxdate): 
                  sys.stdout.write("Oldest sample: %s\n" % seqid)
      output.write(','.join(traits)+'\n')

#

# def writeTipDatesSampling(tipdates, distOutput, opOutput, logOutput, maxDate=None):
#       """
#             Writes distribution and operators for estimating tip dates
#       """

#       for seqid in tipdates.keys():

#             # Max date cannot be in the future
#             if (maxDate != None and maxDate < tipdates[seqid][1]):
#                   tipdates[seqid][1] = maxDate

#             distOutput.write('\n\t\t\t<distribution id="tipDates:%s" monophyletic="false" spec="beast.math.distributions.MRCAPrior" tipsonly="true" tree="@Tree.t:tree">\n' % seqid)
#             distOutput.write('\t\t\t\t<taxonset id="TaxonSet:%s" spec="TaxonSet">\n' % seqid)
#             distOutput.write('\t\t\t\t\t<taxon id="%s" spec="Taxon"/>\n' % seqid)
#             distOutput.write('\t\t\t\t</taxonset>\n')
#             distOutput.write('\t\t\t\t<distr lower="%f" offset="0.0" spec="beast.math.distributions.Uniform" upper="%f"/>\n' % (tipdates[seqid][0], tipdates[seqid][1]))
#             distOutput.write('\t\t\t</distribution>\n')

#             #opOutput.write('\n\t<operator windowSize="1" spec="TipDatesRandomWalker" taxonset="@TaxonSet:%s" tree="@Tree.t:tree" weight="1.0"/>\n' % seqid)
#             #opOutput.write('\n\t<operator windowSize="1" spec="SampledNodeDateRandomWalkerForZeroBranchSATrees" taxonset="@TaxonSet:%s" tree="@Tree.t:tree" weight="1.0"/>\n' % seqid)
#             opOutput.write('\n\t\t<operator windowSize="1" spec="SampledNodeDateRandomWalker" taxonset="@TaxonSet:%s" tree="@Tree.t:tree" weight="1.0"/>\n' % seqid)

#             logOutput.write('\t\t<log idref="tipDates:%s"/>\n' % seqid)
#       #
# #

# def writeDatesCSV(dates, tipdates, output):

#       maxdate = max(dates.values())

#       output.write("id\tdate\tlower\tupper\n")
#       for seqid in dates:

#             if (seqid in tipdates.keys()):
#                   upper = min(tipdates[seqid][1], maxdate)
#                   lower = tipdates[seqid][0]
#             else:
#                   upper = dates[seqid]
#                   lower = dates[seqid]

#             output.write("%s\t%f\t%f\t%f\n" % (seqid, dates[seqid], lower, upper))
#       #
# #


# def getDoubleDate(date): 

#       dec31   = datetime.datetime.strptime(str(date.year)+"-12-31", "%Y-%m-%d")

#       datett  = date.timetuple()
#       dec31tt = dec31.timetuple()

#       return (date.year + float(datett.tm_yday-1)/dec31tt.tm_yday)
# ##



# def getYearDate(datestring, datefmt = "%Y-%m-%d"):
#       """
#       	Takes in date in format 2014-01-01 and returns as year followed by fraction, eg. 2014.0

#    		Also accounts for unknown month or day to return upper and lower bounds
#       """

#       if (datestring.count("-") == 2):
#             date  = datetime.datetime.strptime(datestring, datefmt)
#             upper = date 
#             lower = date 
#       # Only month
#       elif (datestring.count("-") == 1):            
#             date  = datetime.datetime.strptime(datestring+"-15", datefmt)
#             upper = datetime.datetime.strptime("%d-%d-01" % (date.year, date.month), datefmt)
#             if (date.month == 12):
#                   day = 31
#             else:
#                   day = (datetime.date(date.year, date.month+1, 1)-datetime.date(date.year, date.month, 1)).days
#             lower = datetime.datetime.strptime("%d-%d-%d" % (date.year, date.month, day), datefmt)
#       # Only year
#       elif (datestring.count("-") == 0):            
#             date  = datetime.datetime.strptime(datestring+"-07-01", datefmt)
#             upper = datetime.datetime.strptime(datestring+"-01-01", datefmt)
#             lower = datetime.datetime.strptime(datestring+"-12-31", datefmt)
#       else:
#             sys.stdout.write("Unknown date format - %s\n" % datestring)

#       #return ({'date' : getDoubleDate(date), 'upper' : getDoubleDate(upper), 'lower' : getDoubleDate(lower)})
#       return (getDoubleDate(date), getDoubleDate(upper), getDoubleDate(lower))
# #


def makeXMLFile(pars, template, outputpath=""):

	sys.stdout.write(pars["name"]+"...\n")
	formatpars = dict()
	for par in pars:
	     formatpars['$'+par] = pars[par]
	output = template.format(**formatpars)

	if (outputpath == ""):
		outputpath = pars['outputpath']

	if (not os.path.exists(outputpath)):
		os.mkdir(outputpath)

	outfile = open(outputpath+"/"+pars["name"]+".xml", 'w')
	outfile.write(output)
	outfile.close()
#


def formatPars(pars):

      # Read and replace alignments and dates from file
      for par in pars.keys:

            if (par[:9] == "alignment"):

                  #############
                  # Alignment #
                  #############
                  output_align = io.StringIO()

                  if (par.find("_") > 0): 
                        msa = par[par.find("_")+1:]
                  else:
                        msa = "msa"

                  # Write alignment to output_align and return dates to dictionary
                  dates = writeAlignment(pars["alignment"], output_align, msa_id=msa)
                  
                  pars[par]  = output_align.getvalue()
                  output_align.close()

                  #########
                  # Dates #
                  #########
                  output_dates = io.StringIO()

                  if ("ages_"+msa in pars.keys()):
                        pars["ages_"+msa] = open(pars["ages_"+msa]).read().replace("\n",",\n")
                  else:
                        writeDateTrait(dates, output_dates)
                        pars["ages_"+msa] = output_dates.getvalue()

                  output_dates.close()


                  # Redundant, useful for relaxed clock
                  pars["numberOfBranches"] = len(dates) - 2
#


################################################################################################################################  

for filename in sorted(os.listdir(inputpath)):
      if (fnmatch(filename,config)):

            # Load config file
            configfile = open(inputpath+filename, 'r').read().replace("\t"," ")
            pars 	     = yaml.load(configfile)	

            # Set BEAST specific parameters
            outputpath = os.path.abspath(pars["outputpath"] if options.outputpath == '' else options.outputpath)
            template   = open(os.path.abspath(pars["template"] if options.template == '' else options.template), 'r').read()

            # Set parameters not in the config file
            formatPars(pars)

            # Replace config file and save
            makeXMLFile(pars, template, outputpath=outputpath)

            # Output csv of dates
            # outfile = open(outputpath+"/"+pars["name"]+".dates.csv", 'w')
            # writeDatesCSV(pars['dates']['datetrait'], pars['dates']['tipdates'], outfile)
            # outfile.close()
      #
#

######################################################################
# Import a calibration information file and construct a dictionary
######################################################################

def importCalibrationInformation(calfilename):
    """ Import calibration information file

    Return: 
     1. ID of the bank to refine
     2. Calibration Information Dictionary
    """
    try:
        calfile = open(calfilename, "r")
        filelines = calfile.readlines()
        calfile.close()
    except IOError:
        print "Calibration information file %s is not readable or accessible. " % (calfilename)
        raise NotImplementedError("Calibration information file %s is not readable or accessible. " % (calfilename))

    caldict = {}

    bankid = -1
    for rawline in filelines:
        line = rawline.strip()

        if len(line) == 0:
            # Empty line
            continue
        elif line[0] == "#":
            # Comment line
            continue
        else:
            # Information line
            terms = line.split("=")

            if len(terms) != 2:
                # Type of line not defined
                print "Bad line: %s" % (line)
            else:
                # Well defined
                parname = terms[0].strip()
                valuestr = terms[1].strip()

                if parname.lower() == "bank" and valuestr.lower() != "general":
                    # Starting of a new bank
                    bankid = int(valuestr)
                    caldict[bankid] = {}
	        elif parname.lower() == "bank":
		    # Skip information 
		    continue
                else:
                    # Regular Parameter = Value
                    if bankid < 0:
                        caldict[parname] = valuestr
                    else:
                        caldict[bankid][parname] = valuestr
                    # ENDIFELSE
                # ENDIFELSE
            # ENDIFELSE
        # ENNDIFELSE
    # ENDFOR

    return (int(caldict["WORKING_BANKID"]), caldict)

if __name__=="__main__":
	filename = "/home/wzz/Projects/MantidTests/LeBailFit/FinalExam/Calibration_Information.txt"
	wkbankid, caldict = importCalibrationInformation(filename)
	print wkbankid
	print caldict.keys()
	print caldict[1].keys()
	print caldict[1].values()


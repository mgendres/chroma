#
#  $Id: regres.pl,v 1.19 2006-02-28 21:20:05 edwards Exp $
#
#  This is the top-level script used by chroma/scripts/run_chroma_xmldiff.pl
#

sub regresDirs
{
    # 
    # I'd really like to include the regres.pl scripts recursively down 
    # inside the chroma/tests directory, but I'm having difficulty
    # convincing perl to do it. The problem seems to be that the model
    # I want, namely like the c-preprocessor including files that recursively
    # includes other files (in subdirs) is not how the perl "do" works.
    # So, spell out all the many regression dirs and source them individually.
    #
#    return (
#	    "$test_dir/chroma/hadron/mesonspec/regres.pl"
#    );

     return ( 
	    "$test_dir/chroma/gfix/coulgauge/regres.pl",
	    "$test_dir/chroma/glue/fuzwilp/regres.pl",
	    "$test_dir/chroma/glue/wilslp/regres.pl",
	    "$test_dir/chroma/eig/regres.pl",
	    "$test_dir/chroma/hadron/make_source/regres.pl",
	    "$test_dir/chroma/hadron/propagator/regres.pl",
	    "$test_dir/chroma/hadron/sink_smear/regres.pl",
	    "$test_dir/chroma/hadron/seqsource/regres.pl",
	    "$test_dir/chroma/hadron/building_blocks/regres.pl",
	    "$test_dir/chroma/hadron/spectrum/regres.pl",
	    "$test_dir/chroma/hadron/spectrumOct/regres.pl",
	    "$test_dir/chroma/hadron/mesonspec/regres.pl",
	    "$test_dir/chroma/hadron/mres/regres.pl",
	    "$test_dir/chroma/hadron_s/spectrum_s/regres.pl",
	    "$test_dir/chroma/smear/link_smear/regres.pl",
	    "$test_dir/t_leapfrog/regres.pl",
	    "$test_dir/spectrum_s/regres.pl",
	    "$test_dir/hmc/regres.pl",
	    "$test_dir/purgaug/regres.pl"
     );
}

# Return a true value for a require statement
1;

using Pioneer

GetBuildLibParams("/Users/dennisgoldfarb/Downloads/Pioneer Testing/libraries/", 
                    "yeast_prosit",
                    "/Users/dennisgoldfarb/Downloads/Pioneer Testing/FASTA/")
                    #params_path = "/Users/dennisgoldfarb/Downloads/build_test.json")

GetSearchParams("/Users/dennisgoldfarb/Downloads/Pioneer Testing/libraries/yeast_prosit.poin",
                "/Users/dennisgoldfarb/Downloads/Pioneer Testing/MTAC/arrow_out/",
                "/Users/dennisgoldfarb/Downloads/Pioneer Testing/prosit_yeast_results/")
                #params_path = "/Users/dennisgoldfarb/Downloads/search_test.json")

#BuildSpecLib("build/buildspeclib_params.json")
#SearchDIA("build/searchdia_params.json")
#convertMzML("data/mzML/")
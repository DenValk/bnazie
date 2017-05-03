simple_plugin(bnazie "module"
       larcore_Geometry_Geometry_service
       larcore_Geometry
       ${SIMULATIONBASE}
       ${ART_FRAMEWORK_CORE}
       larsim_MCCheater_BackTracker_service
       larsim_MCCheater_PhotonBackTracker_service
       ${ART_FRAMEWORK_PRINCIPAL}
       ${ART_FRAMEWORK_SERVICES_REGISTRY}
       ${ART_FRAMEWORK_SERVICES_OPTIONAL}
       ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
       ${ART_PERSISTENCY_COMMON}
       ${ART_PERSISTENCY_PROVENANCE}
       ${ART_UTILITIES}
       ${MF_MESSAGELOGGER}
       ${MF_UTILITIES}
       ${CETLIB}
       ${ROOT_BASIC_LIB_LIST}

            )
#simple_plugin(m10KTslice "module"
#      Geometry_service
#                        Geometry
#                        BackTracker_service
#                        ${SIMULATIONBASE}
#                        ${ART_FRAMEWORK_CORE}
#      ${ART_FRAMEWORK_PRINCIPAL}
#      ${ART_FRAMEWORK_SERVICES_REGISTRY}
#      ${ART_FRAMEWORK_SERVICES_BASIC}
#                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
#                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
#       ${ART_PERSISTENCY_COMMON}
#      ${ART_PERSISTENCY_PROVENANCE}
#      ${ART_UTILITIES}
#      ${MF_MESSAGELOGGER}
#      ${MF_UTILITIES}
#      ${CETLIB}
#      ${ROOT_BASIC_LIB_LIST}
#              BASENAME_ONLY
#)

install_headers()
install_fhicl()
#install_source( EXTRAS ADDITIONAL_NOTES.txt README.txt )
install_source()
install_scripts()

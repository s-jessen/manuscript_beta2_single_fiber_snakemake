# Rules: renv
#
# All Rules related to the R package `renv`
#
# contributors: @lachlandeer, @julianlanger, @bergmul

# --- renv rules --- #

rule renv_install:
    shell:
        "Rscript -e \"install.packages('renv', repos = 'https://cloud.r-project.org')\""

rule renv_consent:
    shell:
        'Rscript -e "renv::consent(provided = TRUE)"'

rule renv_init:
    shell:
        'Rscript -e "renv::init()"'

rule renv_snap:
    shell:
        'Rscript -e "renv::snapshot()"'

rule renv_restore:
    shell:
        'Rscript -e "renv::restore()"'

rule renv_update:
    shell:
        'Rscript -e "renv::update()"'

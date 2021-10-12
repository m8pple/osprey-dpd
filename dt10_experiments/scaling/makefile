

SIZES:=$(shell seq 12 2 96)
SIZES:=$(foreach s,$(SIZES),$(s)x$(s)x$(s))
MODES=full noangle nobond

DPD_POETS=../../build/dpd-poets


define template

gen/3b16e5m5r0025_$4_$1x$2x$3/dmpci.3b16e5m5r0025_$4_$1x$2x$3 : dmpci.3b16e5m5r0025_$4_WxHxD
	mkdir -p gen/3b16e5m5r0025_$4_$1x$2x$3
	sed -e "s/__W__/$1/g" -e "s/__H__/$2/g" -e "s/__D__/$3/g" $$< > $$@

gen/3b16e5m5r0025_$4_$1x$2x$3/3b16e5m5r0025_$4_$1x$2x$3.1000.state : gen/3b16e5m5r0025_$4_$1x$2x$3/dmpci.3b16e5m5r0025_$4_$1x$2x$3
	( cd gen/3b16e5m5r0025_$4_$1x$2x$3/ && ../../$(DPD_POETS) 3b16e5m5r0025_$4_$1x$2x$3 )

all_initial_states : gen/3b16e5m5r0025_$4_$1x$2x$3/3b16e5m5r0025_$4_$1x$2x$3.1000.state

endef

$(foreach m,$(MODES),$(foreach s,$(SIZES),$(eval $(call template,$(word 1,$(subst x, ,$(s))),$(word 2,$(subst x, ,$(s))),$(word 3,$(subst x, ,$(s))),$(m)))))
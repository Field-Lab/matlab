
cd /mnt/win3/oko/2000-12-12/
y=writeconst('Data000conv',1000,2000,6000,206,65);
y=writeconst('Data001conv',1000,2000,6000,206,65);
y=writeconst('Data002conv',1000,2000,6000,206,65);
y=writeconst('Data003conv',1000,2000,6000,206,65);
y=writeconst('Data006conv',1000,2000,6000,206,65);
y=writeconst('Data008conv',1000,2000,6000,206,65);
y=writeconst('Data009conv',1000,2000,6000,206,65);
y=writeconst('Data010conv',1000,2000,6000,206,65);
y=writeconst('Data011conv',1000,2000,6000,206,65);
y=writeconst('Data012conv',1000,2000,6000,206,65);
y=writeconst('Data013conv',1000,2000,6000,206,65);
y=writeconst('Data014conv',1000,2000,6000,206,65);
%clock

cd /mnt/win3/oko/2000-12-11/
result=convert_data('Data004',206,65,100000);
result=convert_data('Data005',206,65,100000);
%clock
y=writeconst('Data002conv',1000,2000,6000,206,65);
y=writeconst('Data004conv',1000,2000,6000,206,65);
y=writeconst('Data005conv',1000,2000,6000,206,65);
clock
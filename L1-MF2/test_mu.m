gen_data1;

mu = 1;
run_ALM;

mu1_vio_seq = vio_seq;
mu1_obj_seq = obj_seq;
mu1_x = 1:length(mu1_obj_seq);

mu = 3;
run_ALM;

mu3_vio_seq = vio_seq;
mu3_obj_seq = obj_seq;
mu3_x = 1:length(mu3_obj_seq);

mu = 5;
run_ALM;

mu5_vio_seq = vio_seq;
mu5_obj_seq = obj_seq;
mu5_x = 1:length(mu5_obj_seq);

mu = 10;
run_ALM;

mu10_vio_seq = vio_seq;
mu10_obj_seq = obj_seq;
mu10_x = 1:length(mu10_obj_seq);

a = 1;
b = 25;
plot(mu1_x(a:b), mu1_obj_seq(a:b), mu3_x(a:b), mu3_obj_seq(a:b), mu5_x(a:b), mu5_obj_seq(a:b), mu10_x(a:b), mu10_obj_seq(a:b));
legend('\mu = 1', '\mu = 3', '\mu = 5', '\mu = 10');
xlabel('iteration');
ylabel('|M-U\times V''|');

savefig('compare_mu','pdf');

a = 15;
b = 25;
plot(mu1_x(a:b), mu1_obj_seq(a:b), mu3_x(a:b), mu3_obj_seq(a:b), mu5_x(a:b), mu5_obj_seq(a:b), mu10_x(a:b), mu10_obj_seq(a:b));
legend('\mu = 1', '\mu = 3', '\mu = 5', '\mu = 10');
xlabel('iteration');
ylabel('|M-U\times V''|');

savefig('compare_converge_mu','pdf');

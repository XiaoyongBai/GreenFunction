figure ;
hold on;
dis_real=load('./Pak_U_real_CPP.txt');
dis_imag=load('./Pak_U_imag_CPP.txt');

plot(dis_real(:,1), dis_real(:,2));
plot(dis_imag(:,1), dis_imag(:,2));

figure;
hold on;
stress=load('./Pak_SX_CPP.txt');
plot(stress(:,1),stress(:,2));
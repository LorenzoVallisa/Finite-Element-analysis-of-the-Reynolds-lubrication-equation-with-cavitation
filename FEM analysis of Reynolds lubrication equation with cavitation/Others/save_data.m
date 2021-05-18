function save_data(nlevel)

test_name='Test1'
for k=1:nlevel
    [errors,solutions,femregion,Matrices,Dati]=C_main2D(test_name,k);
    name=sprintf('CG_test_level_%i_fem_%s_grid_%s.mat', k, Dati.fem, Dati.type_mesh);
    save(name, 'errors','solutions','femregion','Matrices','Dati');
end
    
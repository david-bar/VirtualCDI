function example1(random_seed, output_path)

rng(random_seed)

random_data = randn(3, 3)

save(output_path, 'random_data')

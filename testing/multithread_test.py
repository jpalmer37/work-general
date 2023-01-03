import multiprocessing as mp
import os
def worker_function(gen, filename):
	"""
	do some work, write results to output
	"""
	#res = f'item: {item} - result: {item ** 2}'

	with open(filename, 'w') as f:
		for i in gen:
			f.write(str(i) + '\n')


if __name__ == '__main__':
	# os.chdir('/home/john.palmer/work/norovirus/')
	print(os.getcwd())

	pool = mp.Pool(4)
	jobs = []
	filenames = ['test.txt' + str(i) for i in range(10)] 

	gens = [[x for x in range(i)] for i in ([10] * 10)]
	
	for gen, file in zip(gens, filenames):
		job = pool.apply_async(worker_function, (gen, file))
		jobs.append(job)

	for job in jobs:
		job.get()

	pool.close()
	pool.join()
from dwave.system import DWaveSampler,EmbeddingComposite
import dwave.inspector

sampler=EmbeddingComposite(DWaveSampler())
h={'s1': 580, 's2': 420, 's3': 144, 's4': 128}
J={('s1','s2'): 152, ('s1','s3'): -144, ('s1','s4'): -512, ('s2','s3'): 16,('s2','s4'):  -512, ('s3','s4'): 128}
sampleset=sampler.sample_ising(h,J,num_reads=100)
for sample in sampleset.samples(sorted_by='energy'):
	print(sample)
dwave.inspector.show(sampleset)

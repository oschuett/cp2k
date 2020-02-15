CFLAGS   = -g -O3 -std=c99 -march=native -Wall -Wextra

all: grid_base_ref_miniapp.x

grid_base_ref_c.o: grid_prepare_pab.o
grid_base_ref_replay.o: grid_base_ref_c.o

%.o: %.c
	$(CC) -c $(CFLAGS) $<

grid_base_ref_miniapp.x: grid_base_ref_miniapp.o grid_base_ref_replay.o grid_base_ref_c.o grid_prepare_pab.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

clean:
	rm -fv grid_base_ref*.o grid_base_ref*.x

#EOF

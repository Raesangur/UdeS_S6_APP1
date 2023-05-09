// GIF642 - Laboratoire - Mémoire partagée inter-processus
// Prépare un espace mémoire partagé et accessible depuis un script Python.
// La synchronisation est effectuée par envoi/réception de messages sur
// stdin/out.
// Messages de diagnostic sur stderr.
// Voir le script associé "waveprop/lab1_ex4.py".
#include <iostream>
#include <sys/shm.h>
#include <sys/mman.h>
#include <sys/fcntl.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <string.h>
#include <array>
#include <vector>
#include <exception>
#include <thread>

// Taille de la matrice de travail (un côté)
static const int MATRIX_SIZE = 100;
static const int BUFFER_SIZE = MATRIX_SIZE * MATRIX_SIZE * MATRIX_SIZE * 3 * sizeof(double);
// Tampon générique à utiliser pour créer le fichier
static char buffer_[BUFFER_SIZE];

using Matrix = std::vector<std::vector<std::vector<std::array<double, 3>>>>;

void wait_signal()
{
    // Attend une entrée (ligne complète avec \n) sur stdin.
    std::string msg;
    std::cin >> msg;
    std::cerr << "CPP: Got signal." << std::endl;
}

void ack_signal()
{
    // Répond avec un message vide.
    std::cout << "" << std::endl;
}

Matrix doH(const Matrix& mtx)
{
    Matrix newMtx{100, std::vector<std::vector<std::array<double, 3>>>{100, std::vector<std::array<double, 3>>{100,}}};
    
    auto thread1 = [](Matrix& newMtx, const Matrix& mtx){
        for(int i = 0; i < MATRIX_SIZE; i++)
        {
            for(int j = 0; j < MATRIX_SIZE - 1; j++)
            {
                for(int k = 0; k < MATRIX_SIZE; k++)
                {
                    newMtx[i][j + 1][k][0] += mtx[i][j + 1][k][2] - mtx[i][j][k][2];
                }
            }
        }

        for(int i = 0; i < MATRIX_SIZE; i++)
        {
            for(int j = 0; j < MATRIX_SIZE; j++)
            {
                for(int k = 0; k < MATRIX_SIZE - 1; k++)
                {
                    newMtx[i][j][k + 1][0] -= mtx[i][j][k + 1][1] - mtx[i][j][k][1];
                }
            }
        }
    };

    auto thread2 = [](Matrix& newMtx, const Matrix& mtx){
        for(int i = 0; i < MATRIX_SIZE; i++)
        {
            for(int j = 0; j < MATRIX_SIZE; j++)
            {
                for(int k = 0; k < MATRIX_SIZE - 1; k++)
                {
                    newMtx[i][j][k + 1][1] += mtx[i][j][k + 1][0] - mtx[i][j][k][0];
                }
            }
        }

        for(int i = 0; i < MATRIX_SIZE - 1; i++)
        {
            for(int j = 0; j < MATRIX_SIZE; j++)
            {
                for(int k = 0; k < MATRIX_SIZE; k++)
                {
                    newMtx[i + 1][j][k][1] -= mtx[i + 1][j][k][2] - mtx[i][j][k][2];
                }
            }
        }
    };


    auto thread3 = [](Matrix& newMtx, const Matrix& mtx){
        for(int i = 0; i < MATRIX_SIZE - 1; i++)
        {
            for(int j = 0; j < MATRIX_SIZE; j++)
            {
                for(int k = 0; k < MATRIX_SIZE; k++)
                {
                    newMtx[i + 1][j][k][2] += mtx[i + 1][j][k][1] - mtx[i][j][k][1];
                }
            }
        }

        for(int i = 0; i < MATRIX_SIZE; i++)
        {
            for(int j = 0; j < MATRIX_SIZE - 1; j++)
            {
                for(int k = 0; k < MATRIX_SIZE; k++)
                {
                    newMtx[i][j + 1][k][2] -= mtx[i][j + 1][k][0] - mtx[i][j][k][0];
                }
            }
        }
    };

    std::thread t1 {thread1, std::ref(newMtx), mtx};
    std::thread t2 {thread2, std::ref(newMtx), mtx};
    thread3(newMtx, mtx);

    t1.join();
    t2.join();

    return newMtx;
}

Matrix doE(const Matrix& mtx)
{
    Matrix newMtx{100, std::vector<std::vector<std::array<double, 3>>>{100, std::vector<std::array<double, 3>>{100,}}};
    
    auto thread1 = [](Matrix& newMtx, const Matrix& mtx){
        for(int i = 0; i < MATRIX_SIZE; i++)
        {
            for(int j = 0; j < MATRIX_SIZE - 1; j++)
            {
                for(int k = 0; k < MATRIX_SIZE; k++)
                {
                    newMtx[i][j][k][0] += mtx[i][j + 1][k][2] - mtx[i][j][k][2];
                }
            }
        }

        for(int i = 0; i < MATRIX_SIZE; i++)
        {
            for(int j = 0; j < MATRIX_SIZE; j++)
            {
                for(int k = 0; k < MATRIX_SIZE - 1; k++)
                {
                    newMtx[i][j][k][0] -= mtx[i][j][k + 1][1] - mtx[i][j][k][1];
                }
            }
        }
    };

    auto thread2 = [](Matrix& newMtx, const Matrix& mtx){
        for(int i = 0; i < MATRIX_SIZE; i++)
        {
            for(int j = 0; j < MATRIX_SIZE; j++)
            {
                for(int k = 0; k < MATRIX_SIZE - 1; k++)
                {
                    newMtx[i][j][k][1] += mtx[i][j][k + 1][0] - mtx[i][j][k][0];
                }
            }
        }

        for(int i = 0; i < MATRIX_SIZE - 1; i++)
        {
            for(int j = 0; j < MATRIX_SIZE; j++)
            {
                for(int k = 0; k < MATRIX_SIZE; k++)
                {
                    newMtx[i][j][k][1] -= mtx[i + 1][j][k][2] - mtx[i][j][k][2];
                }
            }
        }
    };


    auto thread3 = [](Matrix& newMtx, const Matrix& mtx){
        for(int i = 0; i < MATRIX_SIZE - 1; i++)
        {
            for(int j = 0; j < MATRIX_SIZE; j++)
            {
                for(int k = 0; k < MATRIX_SIZE; k++)
                {
                    newMtx[i][j][k][2] += mtx[i + 1][j][k][1] - mtx[i][j][k][1];
                }
            }
        }

        for(int i = 0; i < MATRIX_SIZE; i++)
        {
            for(int j = 0; j < MATRIX_SIZE - 1; j++)
            {
                for(int k = 0; k < MATRIX_SIZE; k++)
                {
                    newMtx[i][j][k][2] -= mtx[i][j + 1][k][0] - mtx[i][j][k][0];
                }
            }
        }
    };

    std::thread t1 {thread1, std::ref(newMtx), mtx};
    std::thread t2 {thread2, std::ref(newMtx), mtx};
    thread3(newMtx, mtx);

    t1.join();
    t2.join();

    return newMtx;
}

void read_matrix(Matrix& mtx, const double* ptr)
{
    const size_t i_off = MATRIX_SIZE * MATRIX_SIZE * 3;
    const size_t j_off = MATRIX_SIZE * 3;
    const size_t k_off = 3;
    
    for(int i = 0; i < MATRIX_SIZE; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            for (int k = 0; k < MATRIX_SIZE; k++)
            {
                mtx[i][j][k][0] = ptr[i * i_off + j * j_off + k * k_off + 0]; 
                mtx[i][j][k][1] = ptr[i * i_off + j * j_off + k * k_off + 1]; 
                mtx[i][j][k][2] = ptr[i * i_off + j * j_off + k * k_off + 2]; 
            }
        }
    }
}

void write_matrix(const Matrix& mtx, double* ptr)
{
    const size_t i_off = MATRIX_SIZE * MATRIX_SIZE * 3;
    const size_t j_off = MATRIX_SIZE * 3;
    const size_t k_off = 3;

    for(int i = 0; i < MATRIX_SIZE; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            for (int k = 0; k < MATRIX_SIZE; k++)
            {
                ptr[i * i_off + j * j_off + k * k_off + 0] = mtx[i][j][k][0];
                ptr[i * i_off + j * j_off + k * k_off + 1] = mtx[i][j][k][1]; 
                ptr[i * i_off + j * j_off + k * k_off + 2] = mtx[i][j][k][2]; 
            }
        }
    }
}


int main(int argc, char** argv)
{
    try
    {
        if (argc < 2)
        {
            std::cerr << "Error : no shared file provided as argv[1]" << std::endl;
            return -1;
        }

        wait_signal();

        // Création d'un fichier "vide" (le fichier doit exister et être d'une
        // taille suffisante avant d'utiliser mmap).
        memset(buffer_, 0, BUFFER_SIZE);
        FILE* shm_f = fopen(argv[1], "w");
        fwrite(buffer_, sizeof(char), BUFFER_SIZE, shm_f);
        fclose(shm_f);

        // On signale que le fichier est prêt.
        std::cerr << "CPP: File ready." << std::endl;
        ack_signal();

        // On ré-ouvre le fichier et le passe à mmap(...). Le fichier peut ensuite
        // être fermé sans problèmes (mmap y a toujours accès, jusqu'à munmap.)
        int shm_fd = open(argv[1], O_RDWR);
        void* shm_mmap = mmap(NULL, BUFFER_SIZE, PROT_WRITE | PROT_READ, MAP_SHARED, shm_fd, 0);
        close(shm_fd);

        if (shm_mmap == MAP_FAILED) {
            std::cerr << "ERROR SHM\n";
            perror(NULL);
            return -1;
        }

        // Pointeur format double qui représente la matrice partagée:
        double* mtx = static_cast<double*> (shm_mmap);
        
        Matrix m{100, std::vector<std::vector<std::array<double, 3>>>{100, std::vector<std::array<double, 3>>{100}}};

        while(true)
        {
            // On attend le signal du parent.
            std::cerr << "CPP: Starting work." << std::endl;

            wait_signal();
            read_matrix(m, mtx);
            Matrix newH = doH(m);
            write_matrix(newH, mtx);
            ack_signal();

            wait_signal();
            read_matrix(m, mtx);
            Matrix newE = doE(m);
            write_matrix(newE, mtx);
            ack_signal();
            
            
            // On signale que le travail est terminé.
            std::cerr << "CPP: Work done." << std::endl;
        }

        munmap(shm_mmap, BUFFER_SIZE);
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << std::endl;
    }
    return 0;
}

################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F90_SRCS += \
../paper2.f90 \
../tprofile.f90 \
../valencia.f90 

OBJS += \
./paper2.o \
./tprofile.o \
./valencia.o 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.f90 subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O0 -g -Wall -c -fmessage-length=0 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

paper2.o: ../paper2.f90

tprofile.o: ../tprofile.f90

valencia.o: ../valencia.f90


clean: clean--2e-

clean--2e-:
	-$(RM) ./paper2.o ./tprofile.o ./valencia.o

.PHONY: clean--2e-


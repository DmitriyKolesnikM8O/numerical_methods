# Компилятор и флаги
CXX := g++
CXXFLAGS := -fsanitize=undefined -pedantic-errors -Wall -Wextra -g3 -std=c++20
INCLUDES = -Itasks/classes
BUILD_DIR := build
TASKS_ROOT := tasks

TASK_DIRS_ROOT := $(wildcard $(TASKS_ROOT)/[0-9]*/)
TASKS_DIRS := $(wildcard $(TASKS_DIRS_ROOT)*[0-9]*.*/)

TASKS := $(shell find $(TASK_DIRS_ROOT) -mindepth 1 -maxdepth 1 -type d -name '*.*' -printf '%P\n')

GENERAL_TASKS := $(filter-out 1.3,$(TASKS))

ALL_TARGETS := $(TASKS) 1.3.2

all: $(ALL_TARGETS)
	@echo "Все задания собраны!"

$(GENERAL_TASKS): %:
	@mkdir -p $(BUILD_DIR)
	@echo "Сборка задания $@..."
	@TASK_ROOT=$$(echo $@ | cut -d. -f1); \
	TASK_PATH=$(TASKS_ROOT)/$$TASK_ROOT/$@; \
	TASK_FILENAME_UNDERSCORE=$$(echo $@ | tr '.' '_'); \
	if [ -f "$$TASK_PATH/$@.cpp" ]; then \
		$(CXX) $(INCLUDES) $(CXXFLAGS) -o $(BUILD_DIR)/$@ $$TASK_PATH/$@.cpp; \
		echo "Собрано: $(BUILD_DIR)/$@ из $$TASK_PATH/$@.cpp"; \
	elif [ -f "$$TASK_PATH/$$TASK_FILENAME_UNDERSCORE.cpp" ]; then \
		$(CXX) $(INCLUDES) $(CXXFLAGS) -o $(BUILD_DIR)/$@ $$TASK_PATH/$$TASK_FILENAME_UNDERSCORE.cpp; \
		echo "Собрано: $(BUILD_DIR)/$@ из $$TASK_PATH/$$TASK_FILENAME_UNDERSCORE.cpp"; \
	else \
		echo "Ошибка: Не найден файл $@.cpp или $$TASK_FILENAME_UNDERSCORE.cpp в $$TASK_PATH"; \
	fi

1.3:
	@mkdir -p $(BUILD_DIR)
	@echo "Сборка задания 1.3 (основной)..."
	@if [ -f "$(TASKS_ROOT)/1/1.3/1_3.cpp" ]; then \
		$(CXX) $(INCLUDES) $(CXXFLAGS) -o $(BUILD_DIR)/1.3 $(TASKS_ROOT)/1/1.3/1_3.cpp; \
		echo "Собрано: $(BUILD_DIR)/1.3 из $(TASKS_ROOT)/1/1.3/1_3.cpp"; \
	else \
		echo "Ошибка: Не найден файл $(TASKS_ROOT)/1/1.3/1_3.cpp"; \
	fi


1.3.2:
	@mkdir -p $(BUILD_DIR)
	@echo "Сборка подзадания 1.3.2..."
	@if [ -f "$(TASKS_ROOT)/1/1.3/1_3_2.cpp" ]; then \
		$(CXX) $(INCLUDES) $(CXXFLAGS) -o $(BUILD_DIR)/1.3.2 $(TASKS_ROOT)/1/1.3/1_3_2.cpp; \
		echo "Собрано: $(BUILD_DIR)/1.3.2 из $(TASKS_ROOT)/1/1.3/1_3_2.cpp"; \
	else \
		echo "Ошибка: Не найден файл $(TASKS_ROOT)/1/1.3/1_3_2.cpp"; \
	fi


list:
	@echo "Доступные задания: $(TASKS)"
	@echo "Доступные подзадания: 1.3.2"
	@echo "Используйте: make <номер_задания>"
	@echo "Например: make 1.5 или make 1.3.2"
	@echo "Или: make all для сборки всего"


clean:
	rm -rf $(BUILD_DIR)

.PHONY: all $(GENERAL_TASKS) 1.3 1.3.2 list clean
